# microscape/cli/metabolism.py
from __future__ import annotations
from pathlib import Path
from typing import Dict, Any, List
import json, csv, typer
from rich.progress import Progress
import numpy as np

from ..io.system_loader import load_system, iter_spot_files_for_env
from ..io.microbe_registry import build_microbe_model_map
from ..io.metabolism_rules import load_rules, MetabolismRules
from ..io.spot_loader import load_spot

app = typer.Typer(add_completion=False, no_args_is_help=True)

def _set_solver(model, solver: str):
    """Best-effort set COBRA solver."""
    try:
        # modern
        model.solver = solver
    except Exception:
        try:
            from cobra.util.solver import change_solver
            change_solver(model, solver)
        except Exception:
            pass

def _apply_bounds_from_conc(model, rules: MetabolismRules, spot_mets: Dict[str, float]) -> List[tuple]:
    """
    For each metabolite in rules.metabolite_map:
      conc = spot_mets.get(Cxxxx, 0)
      lb = -min(conc*mM_per_uptake, max_uptake) if conc>0 else 0
      ub = secretion_upper
    Returns [(rxn_id, lb, ub), ...] actually applied.
    """
    applied = []
    for met_id, ex_id in rules.metabolite_map.items():
        if ex_id not in model.reactions:
            continue
        conc = float(spot_mets.get(met_id, 0.0))
        if rules.uptake.mode != "linear_clip":
            # other modes could be added later
            rate = conc * rules.uptake.mM_per_uptake
            lb = -min(rate, rules.uptake.max_uptake) if conc > 0 else 0.0
        else:
            rate = conc * rules.uptake.mM_per_uptake
            lb = -min(rate, rules.uptake.max_uptake) if conc > 0 else 0.0
        ub = rules.uptake.secretion_upper
        try:
            rxn = model.reactions.get_by_id(ex_id)
            # Ensure lower <= upper
            if lb > ub:
                lb = ub
            rxn.lower_bound = lb
            rxn.upper_bound = ub
            applied.append((ex_id, lb, ub))
        except Exception:
            # skip silently; caller can still solve
            pass
    return applied

@app.command("metabolism")
def metabolism_cmd(
    system_yml: Path = typer.Argument(..., help="Path to system.yml"),
    outdir: Path = typer.Option("outputs/metabolism", help="Output directory"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Extra logs")
):
    """
    Profile steady-state metabolism per spot×microbe:
      - resolve models from system.registry.microbes
      - set exchange bounds from spot metabolite concentrations (mM) via metabolism.yml rules
      - set solver & fallback objective (if model lacks one)
      - run FBA and report objective + selected fluxes
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

    # 2) microbe -> SBML mapping
    microbe_models, warn_models = build_microbe_model_map(sys_info["root"], sys_info["system"].get("registry", {}))
    for w in warn_models:
        typer.secho("WARN: " + w, fg=typer.colors.YELLOW)
    if verbose:
        typer.echo("Resolved models:")
        for k, p in microbe_models.items():
            typer.echo(f"  {k}: {p}")

    # 3) iterate spots and run COBRA
    try:
        import cobra
    except Exception as e:
        typer.secho(f"Cannot import cobra: {e}", fg=typer.colors.RED)
        raise typer.Exit(1)

    env_files = sys_info["environment_files"]
    total_spots = sum(len(list(iter_spot_files_for_env(e, sys_info["paths"]))) for e in env_files) or 1

    all_rows: List[Dict[str, Any]] = []
    per_spot_json: Dict[str, Any] = {}

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

                    # Fallback objective if needed
                    if (model.objective is None or len(model.objective.variables) == 0) and rules.objective:
                        try:
                            if rules.objective in model.reactions:
                                model.objective = model.reactions.get_by_id(rules.objective)
                        except Exception:
                            pass

                    # Per-spot exchange bounds from concentrations
                    applied = _apply_bounds_from_conc(model, rules, spot_mets)

                    # Solve
                    try:
                        sol = model.optimize()
                        status = str(sol.status)
                        obj = float(sol.objective_value) if sol.objective_value is not None else np.nan
                        # record a few exchange fluxes present in the map
                        fx = {}
                        for ex_id in rules.metabolite_map.values():
                            if ex_id in model.reactions and ex_id in sol.fluxes.index:
                                fx[ex_id] = float(sol.fluxes[ex_id])
                        per_microbe[mid] = {"status": status, "objective": obj, "fluxes": fx, "applied": applied}
                    except Exception as e:
                        per_microbe[mid] = {"status": "solve_error", "detail": str(e)}

                    # Flat row for CSV
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

    # 4) write outputs
    (outdir / "metabolism_summary.json").write_text(json.dumps(per_spot_json, indent=2))
    cols = ["spot_id", "microbe", "status", "objective"]
    for r in all_rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    with (outdir / "metabolism_summary.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in all_rows:
            w.writerow(r)

    typer.echo(f"\n🧪 Metabolism profiling complete.")
    typer.echo(f"  Spots processed : {total_spots}")
    typer.echo(f"  Rows written    : {len(all_rows)}")
    typer.echo(f"  CSV             : {(outdir / 'metabolism_summary.csv')}")
    typer.echo(f"  JSON            : {(outdir / 'metabolism_summary.json')}")
