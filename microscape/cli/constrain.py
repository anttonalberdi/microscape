# microscape/cli/constrain.py
from __future__ import annotations
from pathlib import Path
from typing import Dict, Tuple, Optional, List
import json, csv, typer
from rich.progress import Progress
import numpy as np

from ..io.system_loader import load_system, iter_spot_files_for_env, read_spot_yaml, read_microbe_yaml, spot_transcripts, spot_metabolites
from ..io.constraint_rules import load_constraint_rules, env_bounds_from_conc, tx_bounds_from_transcripts, guess_reaction_type
from ..runner.constraints import constrain_one

app = typer.Typer(add_completion=False, no_args_is_help=True)

def _base_bounds_from_model(model) -> Dict[str, Tuple[float,float,str]]:
    out = {}
    for r in model.reactions:
        lb0 = float(r.lower_bound); ub0 = float(r.upper_bound)
        rtype = guess_reaction_type(r)
        out[r.id] = (lb0, ub0, rtype)
    return out

@app.command("constrain")
def constrain_cmd(
    system_yml: Path = typer.Argument(..., help="Path to system.yml"),
    outdir: Path = typer.Option("outputs/constraints", help="Output directory"),
    mode: str = typer.Option("combined", help="environmental|transcriptional|combined"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Extra logs"),
):
    """
    Build constraint proposals per spot×microbe and export:
      - TSV summary (one row per spot×microbe)
      - JSON detailed reactions (one array per spot×microbe)
    """
    sys_info = load_system(system_yml)
    metab_cfg_path = sys_info["configs"]["metabolism"]
    if not metab_cfg_path:
        typer.secho("Metabolism rules not found. Ensure system.config.metabolism points to a file.", fg=typer.colors.RED)
        raise typer.Exit(1)
    rules = load_constraint_rules(metab_cfg_path)

    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    tsv_path = outdir / f"constraints__{mode}.tsv"
    json_path = outdir / f"constraints__{mode}__reactions.json"

    # Build microbe->model path dict
    microbe_models: Dict[str, Path] = {}
    for mid, myml in sys_info["microbe_files"].items():
        d = read_microbe_yaml(mid, sys_info)
        if not d: continue
        model_rel = (d.get("model") or {}).get("path")
        if not model_rel: continue
        mp = (myml.parent / model_rel).resolve()
        if mp.exists():
            microbe_models[mid] = mp
        elif verbose:
            typer.echo(f"WARN: model not found for {mid}: {mp}")

    # iterate
    all_rows: List[Dict] = []
    all_rxn_rows: List[Dict] = []

    env_files = sys_info["environment_files"]
    with Progress() as prog:
        task = prog.add_task("[cyan]Constraining…", total=len(env_files))
        for env_file in env_files:
            for _, spot_path in iter_spot_files_for_env(env_file, sys_info["paths"]):
                spot = read_spot_yaml(spot_path)
                sid = spot.get("name") or spot.get("id") or spot_path.stem
                met_vals = spot_metabolites(spot)
                tr_map   = spot_transcripts(spot)  # {mid:{gene:TPM}}

                microbes_vals = ((spot.get("measurements") or {}).get("microbes") or {}).get("values") or {}
                for mid in microbes_vals.keys():
                    myml = read_microbe_yaml(mid, sys_info)
                    if not myml: 
                        if verbose: typer.echo(f"WARN: Microbe YAML not found for {mid}")
                        continue
                    model_path = (Path(myml["_file"]).parent / (myml.get("model") or {}).get("path","")).resolve()
                    if not model_path.exists():
                        if verbose: typer.echo(f"WARN: Model not found for {mid}: {model_path}")
                        continue

                    import cobra
                    model = cobra.io.read_sbml_model(str(model_path))

                    base_bounds = _base_bounds_from_model(model)
                    env_bounds: Dict[str, Tuple[Optional[float], Optional[float]]] = {}
                    tx_bounds:  Dict[str, Tuple[Optional[float], Optional[float]]] = {}

                    if mode in ("environmental", "combined"):
                        env_bounds = env_bounds_from_conc(met_vals, rules)
                    if mode in ("transcriptional", "combined"):
                        tx_bounds = tx_bounds_from_transcripts(model, tr_map.get(mid, {}), rules)

                    summary_row, rxn_rows = constrain_one(
                        spot_id=sid, microbe_id=mid,
                        base_bounds=base_bounds,
                        env_bounds=env_bounds,
                        tx_bounds=tx_bounds,
                    )
                    summary_row["mode"] = mode
                    all_rows.append(summary_row)
                    all_rxn_rows.extend(rxn_rows)

            prog.advance(task)

    # write TSV summary
    if all_rows:
        cols = ["spot_id","microbe","mode","changed_ex","changed_internal","warnings"]
        with open(tsv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore", delimiter="\t")
            w.writeheader()
            for r in all_rows: w.writerow(r)

    # write JSON full detail
    with open(json_path, "w") as f:
        json.dump(all_rxn_rows, f, indent=2)

    typer.secho("\n🧩 Constraint profiling complete.", fg=typer.colors.GREEN)
    typer.echo(f"  Mode            : {mode}")
    typer.echo(f"  Rows (summary)  : {len(all_rows)}  -> {tsv_path}")
    typer.echo(f"  Reactions (JSON): {len(all_rxn_rows)}  -> {json_path}")
