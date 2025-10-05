# microscape/profile/metabolism.py
from __future__ import annotations
from pathlib import Path
from typing import Dict, Tuple, List, Any

import yaml
from cobra import Model
from cobra.io import read_sbml_model
from cobra.util.solver import linear_reaction_coefficients

from ..io.system_loader import load_system, iter_spot_files_for_env
from ..io.metabolism_rules import MetabolismRules

def _read_yaml(p: Path) -> dict:
    return yaml.safe_load(Path(p).read_text())

def load_microbe_models(system_yml: Path) -> Dict[str, Path]:
    """
    Returns {microbe_id: absolute_path_to_SBML} based on system.registry.microbes entries.
    """
    sys_info = load_system(system_yml)
    root = sys_info["root"]
    registry = (sys_info["system"].get("registry") or {})
    microbes = registry.get("microbes") or []
    out: Dict[str, Path] = {}
    for e in microbes:
        if not isinstance(e, dict):
            continue
        mid = e.get("id")
        myml = root / e.get("file", "")
        if not mid or not myml.exists():
            continue
        mdata = _read_yaml(myml)
        m = mdata.get("microbe") or {}
        model = (m.get("model") or {})
        p = Path(model.get("path", ""))
        sbml = (myml.parent / p) if not p.is_absolute() else p
        if sbml.exists():
            out[mid] = sbml.resolve()
    return out

def _load_spot(spot_yml: Path) -> Tuple[str, Dict[str, float], Dict[str, float]]:
    """
    Returns (spot_id, microbes_counts, metabolite_conc_mM) from a spot YAML.
    """
    d = _read_yaml(spot_yml)
    spot = d.get("spot") or {}
    sid = spot.get("name") or spot.get("id") or spot_yml.stem

    meas = spot.get("measurements") or {}
    microbes = (meas.get("microbes") or {}).get("values") or {}
    microbes = {str(k): float(v) for k, v in microbes.items()}

    mets = (meas.get("metabolites") or {})
    conc = (mets.get("values") or {})  # mM
    conc = {str(k): float(v) for k, v in conc.items()}

    return sid, microbes, conc

def _set_exchange_bounds_from_conc(model: Model, conc_mM: Dict[str, float], rules: MetabolismRules):
    """
    Map spot metabolite concentrations onto exchange reactions in the model.
    Positive uptake is implemented as LOWER bound (negative flux) in COBRA conventions.
    """
    m = rules.metabolite_map  # e.g. {"C0001": "EX_glc__D_e"}
    for met_id, ex_id in m.items():
        if ex_id not in model.reactions:
            continue
        rxn = model.reactions.get_by_id(ex_id)
        c = float(conc_mM.get(met_id, 0.0))
        if c <= 0:
            # no uptake from this pool; allow secretion only
            rxn.lower_bound = 0.0
            rxn.upper_bound = rules.uptake.secretion_upper
        else:
            # convert mM to a simple uptake bound (arbitrary scaling rule)
            uptake = min(c * rules.uptake.mM_per_uptake, rules.uptake.max_uptake)
            # COBRA: uptake is negative flux; allow secretion as well
            rxn.lower_bound = -abs(uptake)
            rxn.upper_bound = rules.uptake.secretion_upper

def _choose_objective(model: Model, objective: str | None):
    """
    Set the objective: 'biomass' (auto-detect), 'biomass_or_sinks', a specific reaction ID,
    or leave model's default if None.
    """
    if not objective:
        return
    obj = str(objective).lower()

    if obj == "biomass":
        # Heuristics: pick reaction with "biomass" in id/name if present
        cand = None
        for r in model.reactions:
            rid = r.id.lower()
            nm = (r.name or "").lower()
            if "biomass" in rid or "biomass" in nm:
                cand = r.id
                break
        if cand:
            model.objective = cand
        return

    if obj == "biomass_or_sinks":
        # Keep current objective if it looks like biomass; else try a known sink exchange
        # If nothing obvious, do nothing (model default)
        coef = linear_reaction_coefficients(model)
        if coef:
            # something already set
            return
        # else try to find any reaction with 'EX_' and positive coeff feasible
        for r in model.reactions:
            if r.id.startswith("EX_"):
                try:
                    model.objective = r.id
                    return
                except Exception:
                    pass
        return

    # treat as explicit reaction id
    if objective in model.reactions:
        model.objective = objective

def _solve_model(model: Model, solver_name: str) -> Tuple[str, float]:
    try:
        model.solver = solver_name
    except Exception:
        # fallback
        pass
    sol = model.optimize()
    status = str(sol.status)
    val = float(getattr(sol, "objective_value", float("nan")))
    return status, val

def profile_spot_metabolism(
    spot_yml: Path,
    microbe_models: Dict[str, Path],
    rules: MetabolismRules,
) -> Tuple[List[dict], dict]:
    """
    For a given spot: run FBA for each microbe with a model, using spot metabolite concentrations
    to set exchange bounds. Returns (rows_for_csv, detail_json_dict).
    """
    sid, microbes_counts, conc_mM = _load_spot(spot_yml)
    rows: List[dict] = []
    details: Dict[str, Any] = {"spot_id": sid, "microbes": {}}

    for mid, count in microbes_counts.items():
        if count <= 0:
            continue
        sbml_path = microbe_models.get(mid)
        if not sbml_path:
            details["microbes"][mid] = {"status": "no_model"}
            continue

        # Load model fresh per microbe/spot (independent media per spot)
        try:
            model = read_sbml_model(str(sbml_path))
        except Exception as e:
            details["microbes"][mid] = {"status": f"sbml_error: {e}"}
            continue

        # Set medium/exchange bounds from spot concentrations
        _set_exchange_bounds_from_conc(model, conc_mM, rules)
        _choose_objective(model, rules.objective)

        status, obj = _solve_model(model, rules.solver)
        rows.append({
            "spot_id": sid,
            "microbe_id": mid,
            "abundance": count,
            "status": status,
            "objective": rules.objective or "model_default",
            "objective_value": obj,
        })
        details["microbes"][mid] = {
            "abundance": count,
            "status": status,
            "objective": rules.objective or "model_default",
            "objective_value": obj,
        }

    return rows, details
