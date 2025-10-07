# microscape/cli/metabolism.py
from __future__ import annotations
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional, Set
import json as jsonlib, csv, typer
from rich.progress import Progress
import numpy as np

from ..io.system_loader import load_system, iter_spot_files_for_env
from ..io.microbe_registry import build_microbe_model_map
from ..io.metabolism_rules import load_rules, MetabolismRules
from ..io.spot_loader import load_spot

app = typer.Typer(add_completion=False, no_args_is_help=True)

def _set_solver(model, solver: str):
    try:
        model.solver = solver
    except Exception:
        try:
            from cobra.util.solver import change_solver
            change_solver(model, solver)
        except Exception:
            pass

def _apply_bounds_from_conc(model, rules: MetabolismRules, spot_mets: Dict[str, float]) -> List[tuple]:
    applied = []
    for met_id, ex_id in rules.metabolite_map.items():
        if ex_id not in model.reactions:
            continue
        conc = float(spot_mets.get(met_id, 0.0))
        rate = conc * rules.uptake.mM_per_uptake
        lb = -min(rate, rules.uptake.max_uptake) if conc > 0 else 0.0
        ub = rules.uptake.secretion_upper
        try:
            rxn = model.reactions.get_by_id(ex_id)
            if lb > ub:
                lb = ub
            rxn.lower_bound = lb
            rxn.upper_bound = ub
            applied.append((ex_id, lb, ub))
        except Exception:
            pass
    return applied

# ---------- constraints support ----------
def _build_constraints_index(detail: dict) -> dict:
    """Normalize a constraints JSON (compact or debug) into:
       idx[spot_id][microbe_id] = {rxn_id: (lb, ub), ...}
    """
    idx: Dict[str, Dict[str, Dict[str, Tuple[Optional[float], Optional[float]]]]] = {}
    if not isinstance(detail, dict):
        return idx
    spots_node = detail.get("spots") or {}

    # Compact shape: spots -> spot_id -> microbes -> mid -> reactions -> rid: {lb,ub}
    compact = bool(spots_node) and all(isinstance(v, dict) and "microbes" in v for v in spots_node.values())
    if compact:
        for sid, s_node in spots_node.items():
            microbes = (s_node or {}).get("microbes") or {}
            for mid, m_node in microbes.items():
                rxns = (m_node or {}).get("reactions") or {}
                for rid, b in rxns.items():
                    if isinstance(b, dict):
                        lb = b.get("lb", b.get("lb_final"))
                        ub = b.get("ub", b.get("ub_final"))
                        if sid not in idx: idx[sid] = {}
                        if mid not in idx[sid]: idx[sid][mid] = {}
                        if lb is not None or ub is not None:
                            try:
                                lbv = float(lb) if lb is not None else None
                                ubv = float(ub) if ub is not None else None
                            except Exception:
                                lbv, ubv = lb, ub
                            idx[sid][mid][rid] = (lbv, ubv)
        return idx

    # Debug shape: spots -> env_id -> spots -> spot_id -> microbes -> mid -> reactions -> rid: {...}
    for _env, env_node in spots_node.items():
        s_spots = (env_node or {}).get("spots") or {}
        for sid, s_node in s_spots.items():
            microbes = (s_node or {}).get("microbes") or {}
            for mid, m_node in microbes.items():
                rxns = (m_node or {}).get("reactions") or {}
                for rid, b in rxns.items():
                    if not isinstance(b, dict):
                        continue
                    lb = b.get("lb", b.get("lb_final"))
                    ub = b.get("ub", b.get("ub_final"))
                    if sid not in idx: idx[sid] = {}
                    if mid not in idx[sid]: idx[sid][mid] = {}
                    if lb is not None or ub is not None:
                        try:
                            lbv = float(lb) if lb is not None else None
                            ubv = float(ub) if ub is not None else None
                        except Exception:
                            lbv, ubv = lb, ub
                        idx[sid][mid][rid] = (lbv, ubv)
    return idx

def _apply_constraints_to_model(model, bounds_map: Dict[str, Tuple[Optional[float], Optional[float]]]) -> Tuple[List[tuple], List[str]]:
    """Apply final bounds to the model. Returns (applied, missing_ids)."""
    applied: List[tuple] = []
    missing: List[str] = []
    if not bounds_map:
        return applied, missing
    for rid, (lb, ub) in bounds_map.items():
        if rid not in model.reactions:
            missing.append(rid)
            continue
        try:
            rxn = model.reactions.get_by_id(rid)
            if lb is not None:
                rxn.lower_bound = float(lb)
            if ub is not None:
                rxn.upper_bound = float(ub)
            applied.append((rid, rxn.lower_bound, rxn.upper_bound))
        except Exception:
            missing.append(rid)
    return applied, missing

def _plan_spots_and_microbes(sys_info) -> Tuple[Set[str], Dict[str, Set[str]]]:
    """Enumerate spots and their microbes for this run."""
    spot_ids: Set[str] = set()
    microbes_by_spot: Dict[str, Set[str]] = {}
    for env_file in sys_info["environment_files"]:
        for sid, spot_path in iter_spot_files_for_env(env_file, sys_info["paths"]):
            spot = load_spot(spot_path)
            spot_id = spot.get("name") or spot.get("id") or sid
            spot_ids.add(spot_id)
            mvals = ((spot.get("measurements") or {}).get("microbes") or {}).get("values") or {}
            microbes_by_spot[spot_id] = set(mvals.keys())
    return spot_ids, microbes_by_spot

@app.command("metabolism")
def metabolism_cmd(
    system_yml: Path = typer.Argument(..., help="Path to system.yml"),
    outdir: Path = typer.Option("outputs/metabolism", help="Output directory"),
    constraints: Path = typer.Option(None, "--constraints", help="constraints__*__reactions.json (compact or debug)"),
    constraints_strict: bool = typer.Option(True, "--constraints-strict/--no-constraints-strict",
                                            help="Error if constraints are unreadable or do not apply to this run."),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Extra logs"),
):
    """
    Profile steady-state metabolism per spot×microbe with optional constraints.
    Writes:
      - legacy: metabolism_summary.csv & metabolism_summary.json (unchanged, no extra echo)
      - named:  metabolism_unconstrained.csv/json OR metabolism_constrained_<mode>.csv/json
    """
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # 1) load system + metabolism rules
    sys_info = load_system(system_yml)
    metab_cfg_path = sys_info.get("metabolism_cfg")
    if not metab_cfg_path:
        typer.secho("Metabolism rules not found. Ensure system.config.metabolism is set.", fg=typer.colors.RED)
        raise typer.Exit(1)
    rules: MetabolismRules = load_rules(metab_cfg_path)

    # Pre-plan the run (used for validating constraints)
    planned_spots, planned_microbes = _plan_spots_and_microbes(sys_info)

    # 2) read & validate constraints (STRICT by default)
    constraints_idx: Dict[str, Dict[str, Dict[str, Tuple[Optional[float], Optional[float]]]]] = {}
    constraints_mode_label: Optional[str] = None  # 'environmental' | 'transcriptional' | 'combined' | 'custom'
    if constraints:
        cpath = Path(constraints)
        if not cpath.exists():
            typer.secho(f"Constraints file not found: {cpath}", fg=typer.colors.RED)
            raise typer.Exit(2)
        try:
            detail = jsonlib.loads(cpath.read_text())
        except Exception as e:
            typer.secho(f"Constraints file is not valid JSON: {e}", fg=typer.colors.RED)
            raise typer.Exit(2)

        # detect mode label for output naming
        raw_mode = str(detail.get("mode", "")).lower().strip()
        if raw_mode in {"environmental", "transcriptional", "combined"}:
            constraints_mode_label = raw_mode
        elif raw_mode:
            constraints_mode_label = "custom"

        constraints_idx = _build_constraints_index(detail)

        if constraints_strict and not constraints_idx:
            typer.secho("Parsed constraints, but found no spot×microbe reactions.", fg=typer.colors.RED)
            raise typer.Exit(2)

        # Check applicability to this run
        applicable_pairs = 0
        applicable_rxns = 0
        for sid, mids in constraints_idx.items():
            if sid not in planned_spots:
                continue
            for mid, rxns in mids.items():
                if mid in planned_microbes.get(sid, set()):
                    applicable_pairs += 1
                    applicable_rxns += len(rxns or {})
        if constraints_strict and (applicable_pairs == 0 or applicable_rxns == 0):
            typer.secho(
                "Constraints loaded but do not apply to any spot×microbe in this run. "
                "Check spot IDs and microbe IDs.", fg=typer.colors.RED
            )
            if verbose:
                typer.echo(f"Planned spots: {sorted(planned_spots)}")
                typer.echo(f"Constraints spots: {sorted(constraints_idx.keys())}")
            raise typer.Exit(2)

        if verbose:
            n_spots = len(constraints_idx)
            n_pairs = sum(len(v) for v in constraints_idx.values())
            typer.echo(f"Loaded constraints for {n_pairs} microbe pairs across {n_spots} constraint spots")

    # 3) microbe -> SBML mapping
    microbe_models, warn_models = build_microbe_model_map(
        sys_info["root"], sys_info["system"], sys_info["paths"]
    )
    for w in warn_models:
        typer.secho("WARN: " + w, fg=typer.colors.YELLOW)
    if verbose:
        typer.echo("Resolved models:")
        for k, p in microbe_models.items():
            typer.echo(f"  {k}: {p}")

    # 4) iterate spots and run COBRA
    try:
        import cobra
    except Exception as e:
        typer.secho(f"Cannot import cobra: {e}", fg=typer.colors.RED)
        raise typer.Exit(1)

    env_files = sys_info["environment_files"]
    total_spots = sum(len(list(iter_spot_files_for_env(e, sys_info["paths"]))) for e in env_files) or 1

    all_rows: List[Dict[str, Any]] = []
    per_spot_json: Dict[str, Any] = {}
    applied_any_constraints = False  # track if at least one bound was set from constraints

    with Progress() as prog:
        task = prog.add_task("[cyan]Profiling metabolism…", total=total_spots)
        for env_file in env_files:
            for sid, spot_path in iter_spot_files_for_env(env_file, sys_info["paths"]):
                spot = load_spot(spot_path)
                spot_id = spot.get("name") or spot.get("id") or sid

                meas = spot.get("measurements", {}) or {}
                spot_mets = (meas.get("metabolites", {}) or {}).get("values", {}) or {}
                spot_microbes = (meas.get("microbes", {}) or {}).get("values", {}) or {}

                per_microbe: Dict[str, Any] = {}

                for mid, _ab in spot_microbes.items():
                    sbml_path = microbe_models.get(mid)
                    if not sbml_path:
                        per_microbe[mid] = {"status": "no_model"}
                        continue
                    # Load SBML
                    try:
                        model = cobra.io.read_sbml_model(str(sbml_path))
                    except Exception as e:
                        per_microbe[mid] = {"status": "sbml_load_error", "detail": str(e)}
                        continue

                    # Solver
                    _set_solver(model, rules.solver)

                    # Fallback objective
                    if (model.objective is None or len(model.objective.variables) == 0) and rules.objective:
                        try:
                            if rules.objective in model.reactions:
                                model.objective = model.reactions.get_by_id(rules.objective)
                        except Exception:
                            pass

                    # Per-spot exchange bounds from concentrations
                    applied_env = _apply_bounds_from_conc(model, rules, spot_mets)

                    # Constraints-based final bounds (override)
                    applied_constraints = []
                    missing_constraints = []
                    if constraints_idx and (spot_id in constraints_idx) and (mid in constraints_idx[spot_id]):
                        applied_constraints, missing_constraints = _apply_constraints_to_model(
                            model, constraints_idx[spot_id][mid]
                        )
                        if applied_constraints:
                            applied_any_constraints = True
                        if constraints_strict and missing_constraints:
                            typer.secho(
                                f"Constraints reference unknown reactions for {spot_id} × {mid}: "
                                f"{', '.join(missing_constraints[:5])}"
                                f"{' …' if len(missing_constraints) > 5 else ''}",
                                fg=typer.colors.RED,
                            )
                            raise typer.Exit(2)

                    # Solve
                    try:
                        with model as mctx:
                            sol = mctx.optimize()
                            status = str(sol.status)
                            obj = float(sol.objective_value) if sol.objective_value is not None else np.nan
                            fx = {}
                            for ex_id in rules.metabolite_map.values():
                                if ex_id in mctx.reactions and ex_id in sol.fluxes.index:
                                    fx[ex_id] = float(sol.fluxes[ex_id])
                            per_microbe[mid] = {
                                "status": status,
                                "objective": obj,
                                "fluxes": fx,
                                "applied_env": applied_env,
                                "applied_constraints": applied_constraints,
                            }
                    except Exception as e:
                        per_microbe[mid] = {"status": "solve_error", "detail": str(e)}

                    # Flat row for CSVs
                    last = per_microbe[mid]
                    row = {
                        "spot_id": spot_id,
                        "microbe": mid,
                        "status": last.get("status"),
                        "objective": last.get("objective"),
                    }
                    for ex_id, val in (last.get("fluxes") or {}).items():
                        row[f"flux:{ex_id}"] = val
                    all_rows.append(row)

                per_spot_json[spot_id] = {"spot_id": spot_id, "microbes": per_microbe}
                prog.advance(task)

    # 5) write outputs
    # named outputs (CSV + JSON), based on constraints type
    if applied_any_constraints:
        label = (constraints_mode_label or "custom")
        base = f"metabolism_constrained_{label}"
    else:
        base = "metabolism_unconstrained"

    named_csv = outdir / f"{base}.csv"
    named_json = outdir / f"{base}.json"

    # write named JSON (same structure as legacy JSON)
    named_json.write_text(jsonlib.dumps(per_spot_json, indent=2))

    # write named CSV (same columns as legacy CSV)
    with named_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in all_rows:
            w.writerow(r)

    # concise report
    typer.echo(f"\n🧪 Metabolism profiling complete.")
    typer.echo(f"  Rows written : {len(all_rows)}")
    typer.echo(f"  Named CSV    : {named_csv}")
    typer.echo(f"  Named JSON   : {named_json}")
