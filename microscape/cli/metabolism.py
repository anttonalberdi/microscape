# microscape/cli/metabolism.py
from __future__ import annotations
from pathlib import Path
from typing import Dict, Any, List
import json, typer
from rich.progress import Progress
import numpy as np

from ..io.system_loader import load_system, iter_spot_files_for_env
from ..io.microbe_registry import build_microbe_model_map
from ..io.metabolism_rules import load_rules, MetabolismRules
from ..io.spot_loader import load_spot  # small helper to read a spot YAML into a dict (see note below)

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command("metabolism")
def metabolism_cmd(
    system_yml: Path = typer.Argument(..., help="Path to system.yml"),
    outdir: Path = typer.Option("outputs/metabolism", help="Output directory"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Extra logs")
):
    """
    Profile steady-state metabolism per spot and microbe using COBRA.
    """
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # 1) Load system + find metabolism rules
    sys_info = load_system(system_yml)
    metab_cfg_path = sys_info.get("metabolism_cfg")
    if not metab_cfg_path:
        typer.secho("Metabolism rules not found. Ensure system.config.metabolism points to a file (relative to config_dir).", fg=typer.colors.RED)
        raise typer.Exit(1)

    rules: MetabolismRules = load_rules(metab_cfg_path)

    # 2) Build microbe -> SBML map
    microbe_models, warn_models = build_microbe_model_map(sys_info["root"], sys_info["system"].get("registry", {}))
    if verbose:
        typer.echo("Resolved microbe models:")
        for mid, p in microbe_models.items():
            typer.echo(f"  - {mid}: {p}")
    for w in warn_models:
        typer.secho("WARN: " + w, fg=typer.colors.YELLOW)

    # 3) Iterate spots and run FBA per microbe present
    env_files = sys_info["environment_files"]
    all_rows: List[Dict[str, Any]] = []
    json_per_spot: Dict[str, Any] = {}

    # Defer heavy imports until needed
    try:
        import cobra
    except Exception as e:
        typer.secho(f"Failed to import cobra: {e}", fg=typer.colors.RED)
        raise typer.Exit(1)

    # Small utility: build default exchange bounds from rules
    # These are global bounds (can be refined with spot metabolite concentrations if desired)
    def default_ex_bounds(met_id: str) -> tuple[float, float]:
        b = rules.uptake.bounds.get(met_id)
        if b is not None:
            return float(b[0]), float(b[1])
        return float(rules.uptake.default_bounds[0]), float(rules.uptake.default_bounds[1])

    # Progress
    total_spots = 0
    for env_file in env_files:
        total_spots += len(iter_spot_files_for_env(env_file, sys_info["paths"]))
    with Progress() as prog:
        task = prog.add_task("[cyan]Profiling metabolism…", total=max(total_spots, 1))

        for env_file in env_files:
            for sid, spot_path in iter_spot_files_for_env(env_file, sys_info["paths"]):
                spot = load_spot(spot_path)  # dict with keys 'name'/'id', 'measurements', etc.
                spot_id = spot.get("name") or spot.get("id") or sid

                # Get metabolite pool at this spot (optional use to adjust bounds later)
                meas = (spot.get("measurements") or {})
                mets_blk = (meas.get("metabolites") or {})
                spot_mets = (mets_blk.get("values") or {})  # {C0001: conc, ...}

                # Microbe abundances (which microbes to consider)
                mic_blk = (meas.get("microbes") or {})
                microbe_counts = (mic_blk.get("values") or {})  # {Mxxxx: count/...}

                per_microbe: Dict[str, Any] = {}

                for mid, abundance in microbe_counts.items():
                    # Skip microbes not in registry or without model
                    sbml_path = microbe_models.get(mid)
                    if sbml_path is None:
                        per_microbe[mid] = {"status": "no_model"}
                        continue

                    # Load model
                    try:
                        model = cobra.io.read_sbml_model(str(sbml_path))
                    except Exception as e:
                        per_microbe[mid] = {"status": "sbml_load_error", "detail": str(e)}
                        continue

                    # Prepare: ensure an objective exists; if not, add a dummy maximize biomass-like sink
                    if model.objective is None or len(model.objective.variables) == 0:
                        # try to set first reaction as objective to avoid errors
                        try:
                            model.objective = list(model.reactions)[-1]
                        except Exception:
                            pass

                    # 3a) Set exchange bounds according to rules (+ optionally scale by spot metabolite availability)
                    # Assumes your exchange reaction IDs are like EX_<met>_e or EX_<met> (rules.ex_mapping maps C0001 -> exchange reaction ID)
                    applied = []
                    for met_id, ex_rxn_id in rules.ex_mapping.items():
                        if ex_rxn_id not in model.reactions:
                            continue
                        rxn = model.reactions.get_by_id(ex_rxn_id)
                        lb, ub = default_ex_bounds(met_id)

                        # If you want to tighten uptake based on spot concentration, do it here (optional):
                        # conc = float(spot_mets.get(met_id, 0.0))
                        # if conc <= rules.uptake.zero_if_missing_below:
                        #     lb = min(lb, 0.0)  # optionally zero
                        # else:
                        #     lb = lb  # keep default or scale

                        try:
                            # cobra uses lower_bound<=upper_bound
                            rxn.lower_bound = lb
                            rxn.upper_bound = ub
                            applied.append((ex_rxn_id, lb, ub))
                        except Exception as e:
                            per_microbe[mid] = {"status": "bound_error", "reaction": ex_rxn_id, "detail": str(e)}
                            break
                    else:
                        # 3b) Solve
                        try:
                            sol = model.optimize()
                            status = str(sol.status)
                            obj = float(sol.objective_value) if sol.objective_value is not None else np.nan

                            # Extract a few fluxes of interest (rules.report_fluxes can list reaction IDs to report)
                            report = {}
                            for rid in rules.report_fluxes:
                                if rid in model.reactions and rid in sol.fluxes.index:
                                    report[rid] = float(sol.fluxes[rid])

                            per_microbe[mid] = {
                                "status": status,
                                "objective": obj,
                                "applied_ex_bounds": applied,
                                "fluxes": report,
                            }
                        except Exception as e:
                            per_microbe[mid] = {"status": "solve_error", "detail": str(e)}

                # Save per spot to JSON dict
                json_per_spot[spot_id] = {
                    "spot_id": spot_id,
                    "microbes": per_microbe,
                }

                # Also collect for flat CSV (one row per spot×microbe)
                for mid, res in per_microbe.items():
                    row = {
                        "spot_id": spot_id,
                        "microbe": mid,
                        "status": res.get("status"),
                        "objective": res.get("objective"),
                    }
                    # flatten a few fluxes if present
                    fx = res.get("fluxes") or {}
                    for k, v in fx.items():
                        row[f"flux:{k}"] = v
                    all_rows.append(row)

                prog.advance(task)

    # 4) Write outputs
    # JSON (per spot, detailed)
    (outdir / "metabolism_summary.json").write_text(json.dumps(json_per_spot, indent=2))

    # CSV
    import csv
    csv_path = outdir / "metabolism_summary.csv"
    # determine all columns seen
    cols: List[str] = ["spot_id", "microbe", "status", "objective"]
    for r in all_rows:
        for k in r.keys():
            if k not in cols:
                cols.append(k)
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in all_rows:
            w.writerow(r)

    # 5) On-screen rollup
    n_rows = len(all_rows)
    n_models = sum(1 for _ in microbe_models.items())
    typer.echo(f"\n🧪 Metabolism profiling complete.")
    typer.echo(f"  Spots processed : {total_spots}")
    typer.echo(f"  Models resolved : {n_models}")
    typer.echo(f"  Rows written    : {n_rows}")
    typer.echo(f"  CSV             : {csv_path}")
    typer.echo(f"  JSON            : {outdir / 'metabolism_summary.json'}")
