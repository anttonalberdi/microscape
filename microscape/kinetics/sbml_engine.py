from __future__ import annotations
from typing import Dict, List
import numpy as np
import warnings
import cobra
from cobra.util.solver import linear_reaction_coefficients

# --- public API ---

def load_models(model_paths: Dict[str, str]) -> Dict[str, cobra.Model]:
    """Load one SBML per guild ID; turn off solver chatter."""
    models = {}
    for gid, p in model_paths.items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            m = cobra.io.read_sbml_model(p)
        m.solver = "glpk"
        models[gid] = m
    return models

def voxel_sources_from_sbml(cfg: dict, models: Dict[str, cobra.Model], C: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    """Return S_f[v]: net production (mmol->conc per unit time) per voxel using expression + Monod."""
    nodes = cfg["space"]["nodes"]
    fmap = cfg["metabolite_map"]                    # field -> exchange rxn id
    fba = cfg.get("fba", {})
    cons = cfg.get("constraints", {})
    substrates = set(fba.get("substrates", []))
    sinks = set(fba.get("sinks", []))
    K = cons.get("uptake_k", {})
    Vmax = cons.get("max_uptake", {})
    alpha = float(cons.get("alpha_flux_to_dC", 1e-3))

    # prep output
    field_names = list(fmap.keys())
    N = len(nodes)
    S = {f: np.zeros(N, float) for f in field_names}

    # per voxel, aggregate across guilds
    for v, nd in enumerate(nodes):
        guilds = nd.get("guilds", {})
        expr_all = nd.get("expression", {})
        # collect contributions
        for gid, frac in guilds.items():
            if frac <= 0: continue
            if gid not in models: continue
            m = models[gid]

            # reset to SBML defaults each call (avoid bound drift)
            reset_bounds_to_sbml(m)

            # Expression -> reaction activity (light E-Flux)
            expr = expr_all.get(gid, {})
            apply_expression_caps(m, expr, floor=0.05)

            # Clamp exchanges by local concentrations (Monod)
            for f, ex_id in fmap.items():
                if ex_id not in m.reactions: continue
                rxn = m.reactions.get_by_id(ex_id)
                # sign convention: uptake is negative (LB negative), secretion positive (UB positive)
                if ex_id in substrates:
                    C_loc = float(C[f][v])
                    cap = monod_cap(C_loc, float(Vmax.get(f, 0.0)), float(K.get(f, 1.0)))
                    rxn.lower_bound = max(-cap, rxn.lower_bound)  # ≤ 0
                    rxn.upper_bound = min(0.0, rxn.upper_bound)   # ≤ 0
                elif ex_id in sinks:
                    # allow secretion; optionally cap with Vmax if provided
                    cap = float(Vmax.get(f, 1e3))
                    rxn.lower_bound = max(0.0, rxn.lower_bound)
                    rxn.upper_bound = min(cap, rxn.upper_bound)

            # Objective: biomass + sinks if present
            set_objective(m, sinks if sinks else None)

            sol = m.optimize()
            if not sol or sol.status != "optimal":
                continue

            # collect exchange fluxes
            for f, ex_id in fmap.items():
                if ex_id in m.reactions:
                    flx = sol.fluxes.get(ex_id, 0.0)  # mmol gDW^-1 h^-1
                    S[f][v] += alpha * frac * flx     # scale by abundance and unit conv.

    return S

# --- helpers ---

def reset_bounds_to_sbml(model: cobra.Model):
    """Reset reaction bounds to SBML-declared defaults; requires latest cobra to keep original bounds."""
    for r in model.reactions:
        if r.lower_bound > r.upper_bound:
            r.lower_bound, r.upper_bound = r.upper_bound, r.lower_bound

def monod_cap(C: float, Vmax: float, K: float) -> float:
    if Vmax <= 0.0: return 0.0
    if C <= 0.0: return 0.0
    if K <= 0.0: return Vmax
    return Vmax * (C / (K + C))

def apply_expression_caps(model: cobra.Model, expr: Dict[str, float], floor: float = 0.05):
    """
    Minimal 'E-Flux' style: if gens are linked by reaction ID (same names),
    scale UB/LB by activity in [floor,1]. If no match, leave as-is.
    You can replace with full GPR parsing later.
    """
    if not expr: return
    # normalize to [0,1] via ranks
    vals = np.array(list(expr.values()), float)
    ranks = np.argsort(np.argsort(vals))
    A = {g: max(float(r)/max(len(vals)-1,1), floor) for g, r in zip(expr.keys(), ranks)}
    for rid, act in A.items():
        if rid in model.reactions:
            rxn = model.reactions.get_by_id(rid)
            scale_bounds(rxn, act, floor)

def scale_bounds(rxn, a: float, floor: float):
    lb, ub = rxn.lower_bound, rxn.upper_bound
    if lb >= 0:            # irreversible forward
        rxn.upper_bound = ub * a
    elif ub <= 0:          # irreversible backward
        rxn.lower_bound = lb * a
    else:                  # reversible
        rxn.lower_bound = lb * a
        rxn.upper_bound = ub * a

def set_objective(model: cobra.Model, sink_ids: List[str] | None):
    # biomass reactions are heuristically recognised by name/id
    biomass = [r for r in model.reactions if "biomass" in r.id.lower()]
    coeffs = {}
    for r in biomass:
        coeffs[r] = 1.0
    if sink_ids:
        for sid in sink_ids:
            if sid in model.reactions:
                coeffs[model.reactions.get_by_id(sid)] = 1.0
    if coeffs:
        model.objective = coeffs
    else:
        # fallback: keep model's default objective
        pass
