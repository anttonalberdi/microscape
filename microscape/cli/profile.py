from __future__ import annotations
from pathlib import Path
from typing import List, Dict, Any
import csv, json, typer
from rich.progress import Progress

from ..io.system_loader import load_system, iter_spot_files_for_env
from ..profile.ecology import load_rules, profile_spot

app = typer.Typer(add_completion=False)

@app.command("profile")
def profile_cmd(
    system_yml: Path = typer.Argument(..., help="Path to system.yml"),
    outdir: Path = typer.Option("outputs/profile", help="Output directory"),
):
    """
    Run ecology profiling over all spots defined in system.yml.
    Produces:
      - profile_summary.csv  (row per spot x microbe)
      - profile_summary.json (same content as JSON)
      - one enriched spot YAML per spot (future extension)
    """
    sys_info = load_system(system_yml)
    root = Path(sys_info["root"])
    ecology_cfg = sys_info["ecology_cfg"]
    env_files: List[Path] = sys_info["environment_files"]

    if not ecology_cfg or not Path(ecology_cfg).exists():
        typer.secho("Ecology rules not found. Ensure system.config.ecology points to a valid file.", fg=typer.colors.RED)
        raise typer.Exit(1)

    rules = load_rules(Path(ecology_cfg))
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Collect rows
    all_rows: List[Dict[str, Any]] = []

    with Progress() as prog:
        task = prog.add_task("[cyan]Profiling…", total=len(env_files) if env_files else 1)
        for env_file in env_files:
            # iterate spots for this env
            for _, spot_path in iter_spot_files_for_env(env_file, sys_info["paths"]):
                rows = profile_spot(spot_path, rules)
                all_rows.extend(rows)
            prog.advance(task)

    # Write CSV & JSON
    out_csv = outdir / "profile_summary.csv"
    out_json = outdir / "profile_summary.json"

    # Determine a stable header (union of keys)
    header_keys: List[str] = []
    for r in all_rows:
        for k in r.keys():
            if k not in header_keys:
                header_keys.append(k)
    # Ensure common columns lead
    lead = [c for c in ["spot", "microbe", "abundance"] if c in header_keys]
    tail = [c for c in header_keys if c not in lead]
    header = lead + tail

    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        for r in all_rows:
            w.writerow(r)

    out_json.write_text(json.dumps(all_rows, indent=2))

    typer.secho("✅ Ecology profiling complete.", fg=typer.colors.GREEN)
    typer.echo(f"CSV : {out_csv}")
    typer.echo(f"JSON: {out_json}")
