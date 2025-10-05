# microscape/profile/metabolism.py
from __future__ import annotations
from pathlib import Path
from typing import Dict, Tuple, List, Any
import json
import csv
import logging

import yaml
import cobra
from cobra.io import read_sbml_model
from cobra.util.solver import linear_reaction_coefficients

from ..io.system_loader import load_system, iter_spot_files_for_env

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

def _read_yaml(p: Path) -> dict:
    return yaml.safe_load(p.read_text())

def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def _silence_cobra_warnings():
    logging.getLogger("cobra").setLevel(logging.ERROR)

def _resolve_microbe_models(system_root: Path, system_obj: dict) -> Dict[str, Path]:
    """
    Build {microbe_id -> model_path} by resolving microbe YAMLs (relative to
    system.paths.microbes_dir) and then resolving microbe.model.path relative
    to each microbe YAML.
    """
    paths = (system_obj.get("paths") or {})
    microbes_dir = Path(paths.get("microbes_dir", "microbes"))
    microbes_dir = (system_root / microbes_dir) if not microbes_dir.is_absolute() else microbes_dir

    models: Dict[str, Path] = {}
    reg = (system_obj.get("registry") or {})
    microbes = reg.get("microbes", []) or []
    for item in microbes:
        if not isinstance(item, dict):
            continue
        mid = item.get("id")
        mfile = item.get("file")
        if not mid or not mfile:
            continue
        microbe_yml = (microbes_dir / mfile) if not Path(mfile).is_absolute() else Path(mfile)
        if not microbe_yml.exists():
            continue
        mdata = _read_yaml(microbe_yml)
        if not isinstance(mdata, dict) or "microbe" not in mdata:
            continue
        model_rel = ((mdata["microbe"].get("model") or {}).get("path"))
        if not model_rel:
            continue
        model_path = (microbe_yml.parent / model_rel) if not Path(model_rel).is_absolute() else Path(model_rel)
        if model_path.exists():
            models[mid] = model_path.resolve()
    return models

def _load_cobra_model_cached(mid: str,
                             model_path: Path,
                             cache: Dict[str, cobra.Model],
                             solver: str) -> cobra.Model:
    if mid in cache:
        return cache[mid]
    # Load model
    model = read_sbml_model(str(model_path))
    # Set solver (GLPK default unless user changed)
    try:
        model.solver = solver
    except Exception:
        # Fallback: keep COBRA default if specified solver isn't available
        pass
    cache[mid] = model
    return model

def _has_objective(model: cobra.Model) -> bool:
    try:
        return bool(linear_reaction_coefficients(model))
    except Exception:
        return False

def _set_objective(model: cobra.Model, reaction_id: str) -> None:
    rxn = model.reactions.get_by_id(reaction_id)
    model.objective = rxn

def _find_exchange(model: cobra.Model, ex_id: str) -> cobra.Reaction | None:
    try:
        return model.reactions.get_by_id(ex_id)
    except KeyError:
        return None

def _set_medium_from_spot(model: cobra.Model,
                          spot_mets: Dict[str, float],
                          metabolite_map: Dict[str, str],
                          uptake_rule: dict) -> Dict[str, Tuple[float, float]]:
    """
    Construct a per-spot medium by setting bounds on specified EX_* reactions.
    Returns the bounds dict we applied so callers can record it: {rxn_id: (lb, ub)}.
    Policy:
      - Start by zeroing uptake (lb=0) for all EX_* we touch (no import by default).
      - For metabolites present in the spot, set lb to -min(max_uptake, mM * mM_per_uptake)
        and ub to secretion_upper.
    """
    applied: Dict[str, Tuple[float, float]] = {}
    mM_per_uptake = float(uptake_rule.get("mM_per_uptake", 1.0))
    max_uptake    = float(uptake_rule.get("max_uptake", 10.0))
    secretion_ub  = float(uptake_rule.get("secretion_upper", 1000.0))

    # First zero out imports for any mapped exchange rxn present in model
    for spot_met, ex_id in metabolite_map.items():
        rxn = _find_exchange(model, ex_id)
        if rxn is None:
            continue
        # default export allowed, import disallowed
        rxn.lower_bound = 0.0
        rxn.upper_bound = secretion_ub
        applied[ex_id] = (rxn.lower_bound, rxn.upper_bound)

    # Then enable uptake where measured concentration > 0
    for spot_met, conc in (spot_mets or {}).items():
        ex_id = metabolite_map.get(spot_met)
        if not ex_id:
            continue
        rxn = _find_exchange(model, ex_id)
        if rxn is None:
            continue
        conc = float(conc)
        if conc <= 0:
            continue
        lb = -min(max_uptake, conc * mM_per_uptake)
        rxn.lower_bound = lb
        rxn.upper_bound = secretion_ub
        applied[ex_id] = (rxn.lower_bound, rxn.upper_bound)

    return applied

# ------------------------------------------------------------------------------
# Main runner
# ------------------------------------------------------------------------------

def run_metabolism(system_yml: str | Path,
                   outdir: str | Path,
                   verbose: bool = False) -> Dict[str, Any]:
    """
    Steady-state FBA profiling per spot and per microbe using COBRApy.

    Writes:
      - metabolism_summary.csv
      - metabolism_summary.json
      - metabolism_summary.info.json
    Returns a small summary dict (also written to JSON).
    """
    _silence_cobra_warnings()

    system_yml = Path(system_yml).resolve()
    outdir = Path(outdir).resolve()
    _ensure_dir(outdir)

    sys_info = load_system(system_yml)
    root: Path = sys_info["root"]
    sysobj: dict = sys_info["system"]
    env_files: List[Path] = sys_info["environment_files"]

    # Resolve metabolism config
    # (mirrors ecology resolution: config.metabolism lives under paths.config_dir unless absolute)
    cfg_rel = (sysobj.get("config") or {}).get("metabolism")
    if not cfg_rel:
        raise FileNotFoundError("system.config.metabolism not defined in system.yml")
    config_dir = (root / (sysobj.get("paths", {}).get("config_dir", "config"))).resolve()
    metab_cfg = (config_dir / cfg_rel) if not Path(cfg_rel).is_absolute() else Path(cfg_rel)
    if not metab_cfg.exists():
        raise FileNotFoundError(f"Metabolism config not found: {metab_cfg}")
    cfg = _read_yaml(metab_cfg).get("metabolism", {})
    solver = str(cfg.get("solver", "glpk"))
    objective_id = cfg.get("objective")  # optional
    metabolite_map: Dict[str, str] = cfg.get("metabolite_map", {}) or {}
    uptake_rule: dict = cfg.get("uptake", {}) or {}

    # Microbe → SBML model path map
    microbe_models = _resolve_microbe_models(root, sysobj)

    # Cache loaded COBRA models by microbe id
    model_cache: Dict[str, cobra.Model] = {}

    # Results containers
    rows: List[Dict[str, Any]] = []
    info = {
        "root": str(root),
        "solver": solver,
        "objective_cfg": objective_id,
        "n_env": 0,
        "n_spot": 0,
        "n_microbe": len(microbe_models),
    }

    # Iterate environments and their spots
    for env_file in env_files:
        env = _read_yaml(env_file).get("environment", {})
        env_id = env.get("id") or env_file.stem
        info["n_env"] += 1

        # gather spot YAMLs using loader (respects env.spots or dirs + system paths)
        spot_list = iter_spot_files_for_env(env_file, sysobj.get("paths", {}))
        if not spot_list:
            continue

        for sid, spot_path in spot_list:
            info["n_spot"] += 1
            spot = _read_yaml(spot_path).get("spot", {})

            meas = spot.get("measurements", {}) or {}
            # microbes
            mblock = meas.get("microbes", {}) or {}
            ab_type = mblock.get("type", "counts")
            abund: Dict[str, float] = (mblock.get("values", {}) or {})
            # metabolites
            met_block = meas.get("metabolites", {}) or {}
            spot_mets: Dict[str, float] = met_block.get("values", {}) or {}

            # Optional: if no microbes measured, skip
            if not abund:
                continue

            # For each microbe present at this spot, run FBA
            for mid, amount in abund.items():
                try:
                    amt = float(amount)
                except Exception:
                    amt = 0.0
                if amt <= 0.0:
                    # still record zero-abundance rows for completeness?
                    # skip for now
                    continue
                model_path = microbe_models.get(mid)
                if not model_path:
                    # microbe present in spot but not in registry/models
                    # record a row to make it visible in the table
                    rows.append({
                        "environment": env_id,
                        "spot": sid,
                        "microbe": mid,
                        "abundance": amt,
                        "abundance_type": ab_type,
                        "status": "no_model",
                        "growth": 0.0
                    })
                    continue

                # Prepare model
                model = _load_cobra_model_cached(mid, model_path, model_cache, solver)

                # Build per-spot medium by setting exchange bounds
                applied_bounds = _set_medium_from_spot(model, spot_mets, metabolite_map, uptake_rule)

                # Objective handling
                if objective_id:
                    try:
                        _set_objective(model, objective_id)
                    except Exception:
                        # keep model's own if missing
                        pass
                # Run optimization
                try:
                    sol = model.optimize()
                    growth = float(sol.objective_value) if sol and sol.status == "optimal" else 0.0
                    status = sol.status if sol else "failed"
                except Exception as e:
                    growth = 0.0
                    status = f"error:{e.__class__.__name__}"

                # Collect selected exchange fluxes (those mapped)
                fluxes: Dict[str, float] = {}
                if model.solution is not None and model.solution.status == "optimal":
                    for smet, ex_id in metabolite_map.items():
                        r = _find_exchange(model, ex_id)
                        if r is None:
                            continue
                        try:
                            v = float(model.solution.fluxes.get(r.id, 0.0))
                        except Exception:
                            v = 0.0
                        fluxes[ex_id] = v

                row = {
                    "environment": env_id,
                    "spot": sid,
                    "microbe": mid,
                    "abundance": amt,
                    "abundance_type": ab_type,
                    "status": status,
                    "growth": growth,
                }
                # add exchange flux columns with stable order
                for ex_id in sorted(metabolite_map.values()):
                    row[f"v_{ex_id}"] = fluxes.get(ex_id, 0.0)

                rows.append(row)

    # ---- Output
    # CSV
    csv_path = outdir / "metabolism_summary.csv"
    if rows:
        # stable columns
        base_cols = ["environment", "spot", "microbe", "abundance", "abundance_type", "status", "growth"]
        ex_cols = sorted({c for r in rows for c in r.keys() if c.startswith("v_EX_") or c.startswith("v_")})
        cols = base_cols + ex_cols
        with csv_path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            for r in rows:
                w.writerow({c: r.get(c, "") for c in cols})
    else:
        # still create an empty file with header
        with csv_path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["environment", "spot", "microbe", "abundance", "abundance_type", "status", "growth"])

    # JSON rows
    json_path = outdir / "metabolism_summary.json"
    json_path.write_text(json.dumps(rows, indent=2))

    # Info
    info_path = outdir / "metabolism_summary.info.json"
    info_path.write_text(json.dumps(info, indent=2))

    return {"n_rows": len(rows), **info}
