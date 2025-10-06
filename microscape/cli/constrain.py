# microscape/cli/constrain.py (header only)
from __future__ import annotations
from pathlib import Path
import json, typer
from rich.progress import Progress

from ..io.system_loader import (
    load_system,
    iter_spot_files_for_env,
    read_spot_yaml,
)
from ..io.metabolism_rules import load_rules
from ..runner.constraints import constrain_one

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command("constrain")
def constrain_cmd(
    system_yml: Path = typer.Argument(..., help="Path to system.yml"),
    outdir: Path = typer.Option("outputs/constrain", help="Output directory"),
):
    outdir.mkdir(parents=True, exist_ok=True)
    sys_info = load_system(system_yml)
    metab_cfg_path = sys_info.get("metabolism_cfg")
    if not metab_cfg_path or not Path(metab_cfg_path).exists():
        typer.secho("No metabolism config found (system.config.metabolism).", fg=typer.colors.RED)
        raise typer.Exit(1)
    rules = load_rules(Path(metab_cfg_path))

    all_rows = []
    env_files = sys_info["environment_files"]

    with Progress() as prog:
        task = prog.add_task("[cyan]Building constraints…", total=max(1, len(env_files)))
        for env_file in env_files:
            for sid, spot_path in iter_spot_files_for_env(env_file, sys_info["paths"]):
                spot = read_spot_yaml(spot_path)
                bounds = constrain_one(spot, rules)
                all_rows.append({
                    "spot_id": sid,
                    "bounds": bounds,
                })
            prog.advance(task)

    # Save JSON for inspection
    (outdir / "constraints.json").write_text(json.dumps(all_rows, indent=2))
    typer.secho(f"✅ Wrote {(outdir / 'constraints.json')}", fg=typer.colors.GREEN)
