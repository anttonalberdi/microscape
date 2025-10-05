# microscape/cli/metabolism.py
from __future__ import annotations
from pathlib import Path
from typing import Dict, List
import json
import typer
from rich.progress import Progress
import csv

from ..io.system_loader import load_system, iter_spot_files_for_env
from ..io.metabolism_rules import load_rules
from ..profile.metabolism import profile_spot_metabolism, load_microbe_models

app = typer.Typer(help="Metabolic steady-state profiling (FBA per spot).", add_completion=False)

@app.command("metabolism")
def metabolism_cmd(
    system_yml: Path = typer.Argument(..., help="Path to system.yml"),
    outdir: Path = typer.Option("outputs/metabolism", help="Output directory"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose errors"),
):
    """
    For every Environment → Spot, run steady-state FBA for each microbe present in the spot,
    using metabolite→exchange mapping and uptake rules from config/metabolism.yml.
    Saves a CSV summary and a JSON with per-spot details.
    """
    sys_info = load_system(system_yml)
    metab_cfg_path = sys_info.get("metabolism_cfg")
    if not metab_cfg_path or not Path(metab_cfg_path).exists():
        typer.secho("Metabolism rules not found. Ensure system.config.metabolism points to a file.", fg=typer.colors.RED)
        raise typer.Exit(1)

    rules = load_rules(Path(metab_cfg_path))

    # Load microbe registry → SBML paths once
    try:
        microbe_models: Dict[str, Path] = load_microbe_models(system_yml)
    except Exception as e:
        typer.secho(f"Failed to map microbe models from system.yml: {e}", fg=typer.colors.RED)
        raise typer.Exit(1)

    env_files = sys_info["environment_files"]
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    all_rows: List[dict] = []
    per_spot_json: Dict[str, dict] = {}

    total_spots = 0
    for env_file in env_files:
        for _, _ in iter_spot_files_for_env(env_file, sys_info["paths"]):
            total_spots += 1

    with Progress() as prog:
        task = prog.add_task("[cyan]Profiling metabolism…", total=total_spots or 1)

        for env_file in env_files:
            for spot_id, spot_path in iter_spot_files_for_env(env_file, sys_info["paths"]):
                try:
                    rows, detail = profile_spot_metabolism(spot_path, microbe_models, rules)
                    all_rows.extend(rows)
                    per_spot_json[spot_id] = detail
                except Exception as e:
                    msg = f"[{spot_id}] {spot_path}: {type(e).__name__}: {e}"
                    if verbose:
                        import traceback
                        typer.secho(msg, fg=typer.colors.RED)
                        typer.echo(traceback.format_exc())
                    else:
                        typer.secho(msg, fg=typer.colors.RED)
                finally:
                    prog.advance(task)

    # Save outputs
    csv_path = outdir / "metabolism_summary.csv"
    if all_rows:
        fieldnames = sorted({k for r in all_rows for k in r.keys()})
        with open(csv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(all_rows)
    else:
        # still write a header to show the run happened
        with open(csv_path, "w", newline="") as f:
            f.write("spot_id,environment_id,microbe_id,obj_value,status\n")

    (outdir / "metabolism_summary.json").write_text(
        json.dumps(per_spot_json, indent=2)
    )

    typer.secho(f"✅ Done. Wrote {csv_path}", fg=typer.colors.GREEN)
