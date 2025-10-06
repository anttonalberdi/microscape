# microscape/cli/constrain.py
from __future__ import annotations
from pathlib import Path
from typing import Dict, List
import json
import typer
from rich.progress import Progress

from ..io.system_loader import (
    load_system,
    iter_spot_files_for_env,
    read_spot_yaml,
)
from ..io.metabolism_rules import load_rules, spot_bounds_from_measurements

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command("constrain")
def constrain_cmd(
    system_yml: Path = typer.Argument(..., help="Path to system.yml"),
    outdir: Path = typer.Option(Path("outputs/constrain"), help="Output directory"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Extra logs"),
):
    """
    Build per-spot exchange bounds from metabolite concentrations using metabolism rules.

    What this does
    --------------
    - resolves environments & spots from system.yml
    - loads metabolism rules (system.config.metabolism)
    - converts each spot's metabolite concentrations (mM) into exchange bounds:
        EX_* lower bound (uptake) and upper bound (secretion)
    - writes a tidy CSV and JSON with one row per (spot × exchange reaction)

    This is *independent* of running FBA; use it to inspect the constraints that
    would be applied during steady-state profiling.
    """
    system_yml = system_yml.resolve()
    outdir = Path(outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    sys_info = load_system(system_yml)
    metab_cfg_path = sys_info.get("metabolism_cfg")
    if not metab_cfg_path or not Path(metab_cfg_path).exists():
        typer.secho(
            "No metabolism rules found. Ensure system.config.metabolism points to a valid YAML.",
            fg=typer.colors.RED,
        )
        raise typer.Exit(1)

    rules = load_rules(Path(metab_cfg_path))
    env_files = sys_info["environment_files"]

    if verbose:
        typer.echo("🧭 Constraint profiling")
        typer.echo(f"  system : {system_yml}")
        typer.echo(f"  rules  : {metab_cfg_path}")
        typer.echo(f"  out    : {outdir}")

    rows: List[Dict] = []
    with Progress() as prog:
        task = prog.add_task("[cyan]Building constraints…", total=max(1, len(env_files)))
        for env_file in env_files:
            for sid, spot_path in iter_spot_files_for_env(env_file, sys_info["paths"]):
                spot = read_spot_yaml(spot_path)
                meas = (spot.get("measurements") or {})
                mets = (meas.get("metabolites") or {})

                bounds = spot_bounds_from_measurements(rules, mets)
                # flatten to rows
                for ex_id, (lb, ub) in bounds.items():
                    rows.append({
                        "spot_id": sid,
                        "exchange_id": ex_id,
                        "lb": float(lb),
                        "ub": float(ub),
                    })
            prog.advance(task)

    # Write CSV
    csv_path = outdir / "constraints.csv"
    if rows:
        import csv
        with csv_path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["spot_id", "exchange_id", "lb", "ub"])
            w.writeheader()
            w.writerows(rows)
    else:
        # still create an empty file with header for consistency
        import csv
        with csv_path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["spot_id", "exchange_id", "lb", "ub"])
            w.writeheader()

    # Write JSON
    json_path = outdir / "constraints.json"
    json_path.write_text(json.dumps(rows, indent=2))

    typer.secho("\n🧱 Constraint profiling complete.", fg=typer.colors.GREEN)
    typer.echo(f"  Spots processed : {sum(1 for _ in rows) and len({r['spot_id'] for r in rows})}")
    typer.echo(f"  Rows written    : {len(rows)}")
    typer.echo(f"  CSV             : {csv_path}")
    typer.echo(f"  JSON            : {json_path}")
