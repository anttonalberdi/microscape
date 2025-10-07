# microscape/cli/ecology.py
from __future__ import annotations
from pathlib import Path
from typing import List
import json, typer
from rich.progress import Progress
from ..io.system_loader import load_system, iter_spot_files_for_env, read_spot_yaml
from ..profile.ecology import profile_spot_ecology  # your working ecology logic

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command("ecology")
def ecology_cmd(
    system_yml: Path = typer.Argument(..., help="Path to system.yml"),
    outdir: Path = typer.Option("outputs/ecology", help="Output directory"),
    verbose: bool = typer.Option(False, "--verbose", "-v")
):
    sys_info = load_system(system_yml)
    root = Path(sys_info["root"])
    ecology_cfg = sys_info.get("ecology_cfg")
    env_files: List[Path] = sys_info["environment_files"]

    if not ecology_cfg or not Path(ecology_cfg).exists():
        typer.secho("Ecology rules not found. Ensure system.config.ecology points to a file.", fg=typer.colors.RED)
        raise typer.Exit(1)

    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    summary_rows = []

    typer.echo("🧭 Profiling ecology")
    typer.echo(f"  system : {system_yml.resolve()}")
    typer.echo(f"  rules  : {ecology_cfg}")
    typer.echo(f"  out    : {outdir.resolve()}")

    with Progress() as prog:
        task = prog.add_task("[cyan]Profiling…", total=len(env_files) or 1)
        for env_file in env_files:
            for sid, spath in iter_spot_files_for_env(env_file, sys_info["paths"]):
                spot = read_spot_yaml(spath).get("spot", {})
                rows = profile_spot_ecology(spot, ecology_cfg)
                # add spot id to rows
                for r in rows:
                    r["spot_id"] = sid
                summary_rows.extend(rows)
            prog.advance(task)

    # write outputs
    (outdir / "profile_summary.json").write_text(json.dumps(summary_rows, indent=2))
    # CSV (flat)
    import csv
    if summary_rows:
        cols = sorted({k for r in summary_rows for k in r.keys()})
        with open(outdir / "profile_summary.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(summary_rows)

    typer.secho("✅ Ecology profiling complete.", fg=typer.colors.GREEN)
