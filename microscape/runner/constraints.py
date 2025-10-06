from __future__ import annotations
from typing import Dict, List, Tuple
from pathlib import Path
import math
import cobra
from cobra import Model
from cobra.io import read_sbml_model, write_sbml_model

from ..io.system_loader import load_system, iter_spot_files_for_env, read_spot_yaml, read_microbe_yaml
from ..io.constraint_rules import ConstraintRules
from ..io.metabolism_rules import (
    MetabolismRules,
    expr_to_exchange_bounds,
    spot_bounds_from_measurements,
)

def constrain_one(
    spot: dict,
    rules: MetabolismRules,
) -> Dict[str, Tuple[float, float]]:
    """
    Build exchange bounds for one spot from its metabolite measurements.
    Returns { exchange_rxn_id: (lb, ub), ... }.
    """
    meas = spot.get("measurements") or {}
    mets = meas.get("metabolites") or {}
    return spot_bounds_from_measurements(rules, mets)

def _normalize_expr(expr_map: Dict[str, float], mode: str, constant: float) -> Dict[str, float]:
    if mode == "constant":
        denom = constant if constant > 0 else 1.0
        return {g: (v/denom) for g,v in expr_map.items()}
    vmax = max(expr_map.values()) if expr_map else 1.0
    if vmax <= 0: vmax = 1.0
    return {g: (v/vmax) for g,v in expr_map.items()}

def _reaction_activity_from_gpr(rx: cobra.Reaction, g2a: Dict[str, float]) -> float:
    """
    Simple rule:
      - If GPR is empty: return 1.0 (unconstrained)
      - AND groups -> min
      - OR between genes -> max
    We rely on reaction.gene_reaction_rule already parsed by cobrapy.
    """
    rule = (rx.gene_reaction_rule or "").strip()
    if not rule:
        # If the model used fbc associations but rule didn't parse, fall back to any gene present
        if rx.genes:
            vals = [g2a.get(g.id, 0.0) for g in rx.genes]
            return max(vals) if vals else 1.0
        return 1.0

    # very simple parser: split by 'or' first, then each token can contain 'and'
    # NOTE: This is a simplification; for complex parentheses use sympy parsing later.
    ors = [t.strip() for t in rule.replace("OR","or").replace("AND","and").split("or")]
    alt_vals = []
    for term in ors:
        and_genes = [g.strip() for g in term.split("and")]
        and_vals = [g2a.get(g, 0.0) for g in and_genes if g]
        if not and_vals:
            alt_vals.append(0.0)
        else:
            alt_vals.append(min(and_vals))
    return max(alt_vals) if alt_vals else 1.0

def _scale_bounds_by_activity(rx: cobra.Reaction, activity: float,
                              min_irrev_ub: float, min_rev_mag: float):
    lb, ub = rx.lower_bound, rx.upper_bound
    # Exchange reactions handled elsewhere—still safe to scale transports/internals
    if lb >= 0.0:  # irreversible forward
        new_ub = max(min_irrev_ub, ub * activity)
        rx.upper_bound = new_ub
        # keep lb (often 0); if you want proportional lb too, do it here
    elif ub <= 0.0:  # irreversible reverse
        new_lb = min(-min_irrev_ub, lb * activity) if min_irrev_ub > 0 else (lb * activity)
        rx.lower_bound = new_lb
    else:  # reversible
        # scale magnitude symmetrically
        rx.lower_bound = lb * activity
        rx.upper_bound = ub * activity
        if min_rev_mag > 0:
            if rx.lower_bound < 0:
                rx.lower_bound = min(rx.lower_bound, -min_rev_mag)
            if rx.upper_bound > 0:
                rx.upper_bound = max(rx.upper_bound,  min_rev_mag)

def _apply_exchange_from_spot(model: Model, spot_mets: Dict[str, float], rules: ConstraintRules):
    if not rules.uptake.enabled:
        return
    fmap = rules.uptake.metabolite_map
    for mid, conc in (spot_mets or {}).items():
        ex_id = fmap.get(mid)
        if not ex_id or ex_id not in model.reactions:
            continue
        rx = model.reactions.get_by_id(ex_id)
        # uptake is negative LB, secretion upper as usual
        uptake_mag = min(rules.uptake.max_uptake, max(0.0, conc * rules.uptake.mM_per_uptake))
        rx.lower_bound = -uptake_mag
        rx.upper_bound = rules.uptake.secretion_upper

def constrain_one(model_path: Path,
                  gene_expr: Dict[str, float],  # gene -> TPM for this microbe at this spot
                  spot_mets: Dict[str, float],  # metabolite id -> concentration
                  rules: ConstraintRules,
                  out_model_path: Path | None) -> Tuple[Model, List[dict]]:
    """
    Returns (constrained_model, per-reaction report rows)
    """
    m = read_sbml_model(str(model_path))

    # Normalize gene expression to 0..1
    norm_mode = rules.gpr_activity.normalize.mode
    const = rules.gpr_activity.normalize.constant
    g_norm = _normalize_expr(gene_expr, norm_mode, const)
    # Threshold + clamp
    th = rules.gpr_activity.threshold_TPM
    cap = rules.gpr_activity.cap_activity
    floor = rules.gpr_activity.floor_activity
    g_act = {g: max(floor, min(cap, (0.0 if v*const < th and norm_mode=="constant" else v))) for g,v in g_norm.items()}
    # (If mode != constant, thresholding is on raw TPM before normalization; adapt if you prefer that convention.)

    # Apply exchange bounds from spot metab (optional)
    _apply_exchange_from_spot(m, spot_mets, rules)

    # Scale internal (and transport) reactions via GPR activity
    rows: List[dict] = []
    for rx in m.reactions:
        if rx.id.startswith("EX_"):
            # handled above; still record activity=1.0 for transparency
            rows.append({
                "reaction": rx.id, "kind": "exchange",
                "lb_old": None, "ub_old": None, "lb_new": rx.lower_bound, "ub_new": rx.upper_bound,
                "activity": 1.0
            })
            continue
        lb_old, ub_old = rx.lower_bound, rx.upper_bound
        a = _reaction_activity_from_gpr(rx, g_act)
        _scale_bounds_by_activity(rx, a, rules.bound_scaling.min_irrev_ub, rules.bound_scaling.min_rev_mag)
        rows.append({
            "reaction": rx.id, "kind": "internal",
            "lb_old": lb_old, "ub_old": ub_old, "lb_new": rx.lower_bound, "ub_new": rx.upper_bound,
            "activity": a
        })

    # Optionally write constrained model
    if out_model_path is not None:
        out_model_path.parent.mkdir(parents=True, exist_ok=True)
        write_sbml_model(m, str(out_model_path))

    return m, rows
