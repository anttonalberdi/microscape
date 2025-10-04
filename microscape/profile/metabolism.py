from __future__ import annotations
from typing import Dict, Any, List, Tuple
from pathlib import Path
import warnings
import yaml

def load_rules(p: Path) -> Dict[str, Any]:
    return yaml.safe_load(Path(p).read_text())

def _linear_uptake(conc_mM: float, k: float, vmax: float) -> float:
    return min(max(conc_mM, 0.0) * max(k, 0.0), vmax)

def _monod_uptake(conc_mM: float, Vmax: float, Km: float) -> float:
    c = max(conc_mM, 0.0)
    return (Vmax * c / (Km + c)) if (Km + c) > 0 else 0.0

def _conc_to_uptake(conc_mM: float, cfg: Dict[str, Any]) -> float:
    up = (cfg.get("uptake") or {})
    model = (up.get("model") or "linear").lower()
    if model == "monod":
        Vmax = float(up.get("monod", {}).get("Vmax", 10.0))
        Km = float(up.get("monod", {}).get("Km", 0.5))
        return _monod_uptake(conc_mM, Vmax, Km)
    # default linear
    k = float(up.get("k_linear", 1.0))
    vmax = float(up.get("vmax", 10.0))
    return _linear_uptake(conc_mM, k, vmax)

def _apply_missing_bounds(model, defaults: Dict[str, float]):
    """Ensure reactions have bounds to silence COBRA 'missing bound' warnings."""
    lb_rev = float(defaults.get("lb_rev", -1000.0))
    ub_rev = float(defaults.get("ub_rev", 1000.0))
    lb_irrev = float(defaults.get("lb_irrev", 0.0))
    ub_irrev = float(defaults.get("ub_irrev", 1000.0))
    for rxn in model.reactions:
        # If bounds undefined or crazy, normalize them
        lb, ub = rxn.lower_bound, rxn.upper_bound
        if lb is None and ub is None:
            # heuristic: exchanges & transports often considered reversible in toy models
            rxn.lower_bound, rxn.upper_bound = lb_rev, ub_rev
        else:
            if rxn.lower_bound is None:
                rxn.lower_bound = lb_irrev if (ub is not None and ub >= 0) else lb_rev
            if rxn.upper_bound is None:
                rxn.upper_bound = ub_rev

def _set_objective(model, obj_cfg: Dict[str, Any]):
    mode = (obj_cfg.get("mode") or "biomass_or_sinks").lower()
    sinks = obj_cfg.get("sinks") or []
    # try to find biomass-like reaction
    biomass = None
    for rxn in model.reactions:
        rid = rxn.id.lower()
        if "biomass" in rid or rid.startswith("biomass"):
            biomass = rxn
            break
    if mode == "biomass" and biomass is not None:
        model.objective = biomass
        return
    if mode == "sinks" or (mode == "biomass_or_sinks" and biomass is None):
        # maximize sum of sink exchange fluxes
        to_max = []
        for ex in sinks:
            if ex in model.reactions:
                to_max.append(model.reactions.get_by_id(ex))
        if to_max:
            model.objective = model.problem.Objective(
                sum(r.flux_expression for r in to_max), direction="max"
            )
            return
    # fallback: biomass if present, else keep whatever default
    if biomass is not None:
        model.objective = biomass

def _apply_medium_from_spot(model, spot_meas: Dict[str, Any], cfg: Dict[str, Any]):
    """
    Set exchange reaction bounds based on spot metabolite concentrations.
    Uptake is negative in COBRA convention.
    """
    metab_map = (cfg.get("metabolism", {}).get("metabolite_map") or {})
    upcfg = (cfg.get("metabolism", {}) or {})
    sec_ub = float(upcfg.get("secretion", {}).get("ub", 1000.0))

    # First, relax all mapped exchanges to no uptake (lb=0), allow secretion up to sec_ub
    for mid, ex_id in metab_map.items():
        if ex_id in model.reactions:
            rxn = model.reactions.get_by_id(ex_id)
            rxn.lower_bound = 0.0
            rxn.upper_bound = sec_ub

    # Then set uptake from concentrations
    metabolites = (spot_meas.get("metabolites") or {})
    concs = metabolites.get("values") or {}
    unit = (metabolites.get("unit") or "mM").lower()
    # If your units change, add a conversion here
    for mid, conc in concs.items():
        ex_id = metab_map.get(mid)
        if not ex_id or ex_id not in model.reactions:
            continue
        vmax = _conc_to_uptake(float(conc), upcfg)
        rxn = model.reactions.get_by_id(ex_id)
        rxn.lower_bound = -float(vmax)  # uptake (negative)
        # rxn.upper_bound left at sec_ub (secretion)

def _collect_exchange_fluxes(model, metab_map: Dict[str, str]) -> Dict[str, float]:
    res = {}
    for mid, ex_id in metab_map.items():
        if ex_id in model.reactions:
            res[mid] = float(model.reactions.get_by_id(ex_id).flux)
    return res

def profile_spot_metabolism(
    spot_yml: Path,
    microbe_models: Dict[str, Path],
    rules: Dict[str, Any],
) -> List[Dict[str, Any]]:
    """
    For each microbe present in the spot measurements, run single-species FBA with the spot medium.
    Returns rows: one per (spot, microbe), including growth and exchange fluxes.
    """
    import cobra
    import contextlib

    sd = yaml.safe_load(Path(spot_yml).read_text())["spot"]
    spot_id = sd.get("name") or sd.get("id") or Path(spot_yml).stem
    meas = (sd.get("measurements") or {})
    microbes = (meas.get("microbes") or {})
    mvals = microbes.get("values") or {}  # {microbe_id: abundance}

    cfg = rules  # already loaded
    metab_map = (cfg.get("metabolism", {}).get("metabolite_map") or {})
    obj_cfg = (cfg.get("metabolism", {}).get("objective") or {})
    defaults = (cfg.get("metabolism", {}).get("defaults") or {})

    rows: List[Dict[str, Any]] = []

    for m_id, abundance in mvals.items():
        model_path = microbe_models.get(m_id)
        if not model_path:
            # Skip microbes with no model path
            rows.append({
                "spot": spot_id, "microbe": m_id, "abundance": abundance,
                "status": "no_model", "growth": 0.0
            })
            continue

        # Load model (quiet known warnings)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = cobra.io.read_sbml_model(str(model_path))

        # Normalize missing bounds to avoid spam
        _apply_missing_bounds(model, defaults)
        _set_objective(model, obj_cfg)
        _apply_medium_from_spot(model, meas, cfg)

        # Solve
        sol = model.slim_optimize(error_value=None)  # returns float or None
        status = model.solver.status  # string
        growth = float(sol) if sol is not None else 0.0

        # Get exchange fluxes for mapped metabolites
        ex_flux = _collect_exchange_fluxes(model, metab_map)

        row = {
            "spot": spot_id,
            "microbe": m_id,
            "abundance": abundance,
            "status": status,
            "growth": growth,
        }
        # add metabolite-labeled fluxes (positive = secretion, negative = uptake)
        for mid, v in ex_flux.items():
            row[f"flux_{mid}"] = round(v, 6)
        rows.append(row)

    return rows
