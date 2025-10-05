from __future__ import annotations
from pathlib import Path
from typing import List, Dict, Any
import csv, json, typer
from rich.progress import Progress

from ..io.system_loader import load_system, iter_spot_files_for_env
from ..profile.metabolism import run_metabolism

app = typer.Typer(add_completion=False)

@app.command("metabolism")
def metabolism_cmd(
    system_yml: Path = typer.Argument(..., help="Path to system.yml"),
    outdir: Path = typer.Option("outputs/profile_metabolism", help="Output directory"),
):
    """
    Run metabolic profiling (FBA per microbe per spot) using SBML models and spot metabolite concentrations.
    Outputs:
      - metabolism_summary.csv / .json
    """
    sys_info = load_system(system_yml)
    root = Path(sys_info["root"])
    metab_cfg = (sys_info["system"].get("config") or {}).get("metabolism")
    config_dir = root / (sys_info["paths"].get("config_dir") or "config")
    metab_cfg_path = (Path(metab_cfg) if Path(metab_cfg).is_absolute() else (config_dir / metab_cfg))

    if not metab_cfg_path.exists():
        typer.secho("Metabolism rules not found. Ensure system.config.metabolism points to a valid file.", fg=typer.colors.RED)
        raise typer.Exit(1)

    rules = load_rules(metab_cfg_path)

    # Build microbe -> SBML path map from system registry
    microbe_models: Dict[str, Path] = {}
    for entry in (sys_info["system"].get("registry", {}).get("microbes") or []):
        # entries are {id: M0001, file: M0001.yml}
        myml = root / entry.get("file")
        if not myml.exists():
            continue
        m = (yaml.safe_load(myml.read_text()) or {}).get("microbe", {})
        model_path = (m.get("model") or {}).get("path")
        if model_path:
            p = Path(model_path)
            if not p.is_absolute():
                p = myml.parent / p
            microbe_models[entry["id"]] = p

    env_files: List[Path] = sys_info["environment_files"]
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)

    all_rows: List[Dict[str, Any]] = []
    with Progress() as prog:
        # count spots for progress total
        n_spots = 0
        for env_file in env_files:
            n_spots += len(iter_spot_files_for_env(env_file, sys_info["paths"]))
        task = prog.add_task("[cyan]Profiling metabolism…", total=max(n_spots, 1))

        for env_file in env_files:
            for _, spot_path in iter_spot_files_for_env(env_file, sys_info["paths"]):
                rows = profile_spot_metabolism(spot_path, microbe_models, rules)
                all_rows.extend(rows)
                prog.advance(task)

    # Save
    out_csv = outdir / "metabolism_summary.csv"
    out_json = outdir / "metabolism_summary.json"

    # Make a stable header
    headers = []
    for r in all_rows:
        for k in r.keys():
            if k not in headers:
                headers.append(k)
    lead = [c for c in ["spot", "microbe", "abundance", "status", "growth"] if c in headers]
    tail = [c for c in headers if c not in lead]
    header = lead + tail

    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        for r in all_rows:
            w.writerow(r)

    out_json.write_text(json.dumps(all_rows, indent=2))

    typer.secho("✅ Metabolic profiling complete.", fg=typer.colors.GREEN)
    typer.echo(f"CSV : {out_csv}")
    typer.echo(f"JSON: {out_json}")
