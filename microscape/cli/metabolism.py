# microscape/cli/metabolism.py
from __future__ import annotations
from pathlib import Path
from typing import Dict, List
import json, typer, csv
from rich.progress import Progress
from ..io.system_loader import load_system, iter_spot_files_for_env, read_spot_yaml, read_microbe_yaml
from ..io.metabolism_rules import load_rules, MetabolismRules
from ..runner.metabolism import profile_spot_metabolism  # your steady-state FBA function

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command("metabolism")
def metabolism_cmd(
    system_yml: Path = typer.Argument(..., help="Path to system.yml"),
    outdir: Path = typer.Option("outputs/metabolism", help="Output directory"),
    verbose: bool = typer.Option(False, "--verbose", "-v")
):
    sys_info = load_system(system_yml)
    metab_cfg_path = sys_info.get("metabolism_cfg")
    if not metab_cfg_path or not Path(metab_cfg_path).exists():
        typer.secho("Metabolism rules not found. Ensure system.config.metabolism is set.", fg=typer.colors.RED)
        raise typer.Exit(1)

    rules: MetabolismRules = load_rules(Path(metab_cfg_path))
    microbe_models: Dict[str, Path] = {}
    for mid, myml in sys_info["microbe_registry"].items():
        md = read_microbe_yaml(mid, sys_info)
        if not md or "microbe" not in md:
            if verbose: typer.echo(f"WARN: Microbe YAML not found: {myml}")
            continue
        model_rel = (md["microbe"].get("model") or {}).get("path")
        if not model_rel:
            continue
        mp = (Path(myml).parent / model_rel).resolve()
        if mp.exists():
            microbe_models[mid] = mp

    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    env_files: List[Path] = sys_info["environment_files"]

    typer.echo("Profiling metabolism…")
    with Progress() as prog:
        total = sum(len(iter_spot_files_for_env(e, sys_info["paths"])) for e in env_files) or 1
        task = prog.add_task("[cyan]Profiling metabolism…", total=total)
        for env_file in env_files:
            for _, spot_path in iter_spot_files_for_env(env_file, sys_info["paths"]):
                spot = read_spot_yaml(spot_path).get("spot", {})
                rows.extend(profile_spot_metabolism(spot, microbe_models, rules))
                prog.advance(task)

    # write
    (outdir / "metabolism_summary.json").write_text(json.dumps(rows, indent=2))
    if rows:
        cols = sorted({k for r in rows for k in r.keys()})
        with open(outdir / "metabolism_summary.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)

    typer.echo("\n🧪 Metabolism profiling complete.")
    typer.echo(f"  Spots processed : {len(set(r['spot_id'] for r in rows)) if rows else 0}")
    typer.echo(f"  Rows written    : {len(rows)}")
    typer.echo(f"  CSV             : {outdir / 'metabolism_summary.csv'}")
    typer.echo(f"  JSON            : {outdir / 'metabolism_summary.json'}")
