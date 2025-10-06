from __future__ import annotations
from pathlib import Path
from typing import Dict, List
import csv, json, typer
from rich.progress import Progress

from ..io.system_loader import load_system, iter_spot_files_for_env, read_spot_yaml, read_microbe_yaml
from ..io.constraint_rules import load_constraint_rules, ConstraintRules
from ..runner.constraints import constrain_one

app = typer.Typer(help="Apply gene expression & environment constraints to GSMMs (per spot).")

@app.command()
def constrain(
    system_yml: Path = typer.Argument(..., help="Path to system.yml"),
    outdir: Path = typer.Option("outputs/constraints", help="Where to write summaries/models"),
    write_models: bool = typer.Option(True, help="Write contextualized SBMLs")
):
    sys_info = load_system(system_yml.resolve())
    rules_path = (sys_info["system"].get("config") or {}).get("constraints")
    if not rules_path:
        typer.secho("constraints config not set in system.config.constraints", fg=typer.colors.RED)
        raise typer.Exit(1)

    # Resolve constraints.yml under config_dir
    config_dir = (sys_info["paths"] or {}).get("config_dir")
    cfg_base = sys_info["root"]
    if config_dir:
        cfg_base = (sys_info["root"] / config_dir)
    rules_file = (cfg_base / rules_path) if not Path(rules_path).is_absolute() else Path(rules_path)
    rules: ConstraintRules = load_constraint_rules(rules_file)

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    model_out_base = outdir / (rules.write_models.dir if write_models and rules.write_models.enabled else "models_disabled")

    # Build microbe -> model path map
    microbe_models: Dict[str, Path] = {}
    for mreg in (sys_info["system"].get("registry") or {}).get("microbes", []):
        if not isinstance(mreg, dict): continue
        mfile = (sys_info["root"] / (sys_info["paths"].get("microbes_dir","microbes")) / mreg.get("file")).resolve()
        mdata = read_microbe_yaml(mfile)
        mp = (mfile.parent / (mdata["microbe"]["model"]["path"])).resolve()
        microbe_models[mdata["microbe"]["id"]] = mp

    rows: List[dict] = []
    env_files = sys_info["environment_files"]
    # Count spots
    total_spots = 0
    for env_file in env_files:
        for _sid, _sp in iter_spot_files_for_env(env_file, sys_info["paths"]):
            total_spots += 1

    with Progress() as prog:
        task = prog.add_task("[cyan]Constraining models…", total=total_spots)
        for env_file in env_files:
            env_dir = env_file.parent
            env = read_spot_yaml(env_file)  # reusing helper name; just loads YAML
            env_id = env["environment"]["id"]
            for sid, spot_path in iter_spot_files_for_env(env_file, sys_info["paths"]):
                spot = read_spot_yaml(spot_path)["spot"]
                # transcripts
                tx = (spot.get("measurements") or {}).get("transcripts") or {}
                per_microbe = (tx.get("values") or {})
                # metabolites
                mets = ((spot.get("measurements") or {}).get("metabolites") or {}).get("values", {})

                for mid, model_path in microbe_models.items():
                    gene_expr = per_microbe.get(mid) or {}
                    if not gene_expr:
                        # no expression for this microbe at this spot—skip or set all-zero?
                        continue
                    # output model file path
                    out_model_path = None
                    if write_models and rules.write_models.enabled:
                        out_model_path = (model_out_base / env_id / sid / f"{mid}.xml").resolve()
                    model, rx_rows = constrain_one(model_path, gene_expr, mets, rules, out_model_path)
                    # Append summarized rows (one per reaction) with keys for joining later
                    for r in rx_rows:
                        r.update({"env_id": env_id, "spot_id": sid, "microbe": mid})
                        rows.append(r)
                prog.advance(task, 1)

    # Write summary
    csv_path = outdir / "constraints_summary.csv"
    json_path = outdir / "constraints_summary.json"
    fieldnames = ["env_id","spot_id","microbe","reaction","kind",
                  "lb_old","ub_old","lb_new","ub_new","activity"]
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})
    json_path.write_text(json.dumps({"n_rows": len(rows), "rows": rows[:1000]}, indent=2))
    typer.secho(f"✅ Constrained {len(rows)} reaction-rows across spots. CSV: {csv_path}", fg=typer.colors.GREEN)
