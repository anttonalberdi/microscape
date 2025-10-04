from __future__ import annotations
from pathlib import Path
import json, typer, shutil, sys, subprocess
import numpy as np
from rich.progress import Progress
from ..runner.snapshot import run_snapshot
from ..io.graph_config import load_graph_yaml
from ..viz.graph import scatter_field, interpolate_to_grid
from ..validation.system import validate_system
from ..io.system_loader import load_system, iter_spot_files_for_env
from ..profile.ecology import load_rules, profile_spot

app = typer.Typer(add_completion=False, no_args_is_help=True)

# Repository URL used for updates
REPO_URL = "git+https://github.com/anttonalberdi/microscape.git"

@app.command()
def update(
    ref: str = typer.Option("main", help="Git ref: branch, tag, or commit."),
    yes: bool = typer.Option(False, "--yes", "-y", help="Skip confirmation."),
    with_deps: bool = typer.Option(False, help="Reinstall dependencies too."),
    quiet: bool = typer.Option(True, "--quiet/--no-quiet", help="Reduce pip output.", show_default=True),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Show full pip logs for debugging."),
):
    """
    Update MicroScape from GitHub. Tries a fast path first (no deps, no build isolation).
    If the build backend isn't present (e.g. hatchling), it retries with build isolation.
    """
    url = f"{REPO_URL}@{ref}"
    if not yes:
        typer.confirm(f"This will install from {url} into ({sys.executable}). Continue?", abort=True)

    base = [sys.executable, "-m", "pip", "install", "--upgrade", "--force-reinstall"]
    if quiet and not verbose:
        base.append("-q")

    def run_cmd(cmd, capture=False):
        cmd_str = " ".join(cmd)
        if verbose:
            typer.echo("Running: " + cmd_str)
            return subprocess.run(cmd, check=True)
        if capture:
            # capture output to avoid scary tracebacks in fast path
            return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        else:
            return subprocess.run(cmd, check=True)

    # 0) Nice hint if hatchling is missing (helps fast path succeed)
    if not verbose and shutil.which("python"):
        try:
            __import__("hatchling.build")
        except Exception:
            typer.echo("ℹ️  Tip: install 'hatchling' to enable faster updates (pip install hatchling).")

    # 1) Fast path: no deps, no build isolation (captured to suppress traceback noise)
    fast_cmd = base + ["--no-deps", "--no-build-isolation", url]
    if verbose:
        typer.echo("Trying fast update (no-deps, no-build-isolation)…")
    res = run_cmd(fast_cmd, capture=not verbose)
    if res.returncode == 0:
        typer.secho("MicroScape updated (fast path).", fg=typer.colors.GREEN)
        return

    # If fast path failed but we captured output, keep it for later if needed
    if not verbose:
        typer.secho("Fast update failed; retrying with build isolation…", fg=typer.colors.YELLOW)

    # 2) Retry with build isolation; include deps only if requested
    retry_cmd = base + ([url] if not with_deps else [url])
    try:
        run_cmd(retry_cmd)
        typer.secho("MicroScape updated (build isolation).", fg=typer.colors.GREEN)
    except subprocess.CalledProcessError as e:
        if not verbose:
            typer.secho("❌ Update failed. You can try:", fg=typer.colors.RED)
            typer.echo("   • pip install hatchling")
            typer.echo("   • microscape update --with-deps")
            typer.echo("   • microscape update --verbose   # to see full pip logs")
        raise typer.Exit(code=e.returncode)

"""
Validate a microscape project starting at system.yml.
Exits with non-zero code if errors are found.
"""

@app.command("validate")
def validate_cmd(
    system_yml: Path = typer.Argument(..., help="Path to system.yml"),
    json_out: bool = typer.Option(False, "--json", help="Emit JSON report"),
):

    summary, errors, warnings = validate_system(system_yml.resolve())
    if json_out:
        typer.echo(json.dumps({"summary": summary, "errors": errors, "warnings": warnings}, indent=2))
    else:
        typer.echo("== microscape validation report ==")
        typer.echo(f"Root: {summary.get('root')}")
        typer.echo(f"Microbes: {summary.get('microbes',{}).get('count',0)}")
        typer.echo(f"Environments: {summary.get('environments',{}).get('count',0)}")
        typer.echo(f"Spots: {summary.get('spots',{}).get('count',0)}")
        models = summary.get('models', {})
        space = summary.get('space', {})
        typer.echo(f"Metabolites (unique): {models.get('metabolites_unique', 0)} | Species (total): {models.get('species_total', 0)}")
        typer.echo(f"Genes: total {models.get('genes_total',0)}, expressed {models.get('genes_expressed',0)}, unused {models.get('genes_unused',0)}")
        dims = space.get('dimensions')
        dim_label = 'none' if dims in (None, 0) else (f"{dims}D")
        typer.echo(f"Spatial positions: {'yes' if space.get('has_positions') else 'no'}; spots with pos: {space.get('spots_with_position',0)}; dims: {dim_label}")
        if warnings:
            typer.echo("")
            typer.secho(f"Warnings ({len(warnings)}):", fg=typer.colors.YELLOW)
            for w in warnings:
                typer.echo(f"  - {w}")
        if errors:
            typer.echo("")
            typer.secho(f"Errors ({len(errors)}):", fg=typer.colors.RED)
            for e in errors:
                typer.echo(f"  - {e}")
        typer.echo("")
        typer.secho("Result: " + ("FAILED" if errors else "OK"), fg=typer.colors.RED if errors else typer.colors.GREEN)
    raise typer.Exit(code=1 if errors else 0)

"""
    Run ecology profiling on a microscape project:
    - reads system.yml -> environments -> spots
    - infers per-microbe ecology states/scores from transcripts (TPM) using ecology config
    - writes enriched spot YAMLs under OUTDIR mirroring the original layout
    - writes profile_summary.csv (and .json if requested)
"""

@app.command("profile")
def profile_cmd(
    system_yml: Path = typer.Argument(..., help="Path to system.yml"),
    outdir: Path = typer.Option("outputs/profile", help="Directory to write enriched outputs"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Print discovered files"),
):

    sys_info = load_system(system_yml)
    root = sys_info["root"]
    env_files = sys_info["environment_files"]
    rules_path = sys_info["ecology_cfg"]
    if not rules_path or not rules_path.exists():
        typer.secho("❌ Ecology rules not found (system.config.ecology).", fg=typer.colors.RED)
        raise typer.Exit(2)

    rules = load_rules(rules_path)
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / "profile_summary.csv"
    json_path = outdir / "profile_summary.json"

    typer.echo("🧭 Profiling ecology")
    typer.echo(f"  system : {system_yml.resolve()}")
    typer.echo(f"  rules  : {rules_path}")
    typer.echo(f"  out    : {outdir.resolve()}")

    headers = ["spot","microbe","abundance"]
    # infer trait columns from rules
    for t in (rules.get("ecology",{}).get("traits") or []):
        tid = t["id"]
        headers += [f"{tid}_state", f"{tid}_score", f"{tid}_expr_TPM"]

    total_rows = 0
    env_count = 0
    with Progress() as progress, open(csv_path, "w", newline="") as f:
        task = progress.add_task("[cyan]Profiling…", total=None)
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()

        for env_file in env_files:
            env_count += 1
            spot_list = iter_spot_files_for_env(env_file, sys_info["paths"])
            if verbose:
                typer.echo(f"  • {env_file.name}: {len(spot_list)} spots")
            for sid, spath in spot_list:
                if not spath or not spath.exists():
                    if verbose:
                        typer.echo(f"    - SKIP missing spot file: {spath}")
                    continue
                rows = profile_spot(spath, rules)
                for r in rows:
                    writer.writerow(r)
                    total_rows += 1
                progress.advance(task)

    summary = {
        "n_environments": env_count,
        "n_spots_processed": total_rows,  # rows count (microbe-spot pairs)
        "traits": [t["id"] for t in (rules.get("ecology",{}).get("traits") or [])]
    }
    json_path.write_text(json.dumps(summary, indent=2))
    typer.secho("✅ Ecology profiling complete.", fg=typer.colors.GREEN)


@app.command("simulate")
def simulate_cmd(
    config: Path = typer.Argument(..., help="Graph YAML config"),
    outdir: Path = typer.Option("outputs/snapshot", help="Output directory"),
    plot: bool = typer.Option(True, help="Save plots"),
):
    outdir.mkdir(parents=True, exist_ok=True)
    with Progress() as prog:
        task = prog.add_task("[cyan]Evaluating snapshot…", total=100)
        R, summary = run_snapshot(str(config))
        prog.update(task, completed=100)
    np.savez_compressed(outdir / "residuals.npz", **R)
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    if plot:
        cfg = load_graph_yaml(config)
        nodes = cfg["space"]["nodes"]; edges = cfg["space"]["edges"]
        pos = np.array([nd.get("pos_um", [0,0])[:2] for nd in nodes], float)
        id_to_idx = {nd["id"]: i for i, nd in enumerate(nodes)}
        edge_index = np.array([[id_to_idx[e["i"]], id_to_idx[e["j"]]] for e in edges], int)
        for k, arr in R.items():
            scatter_field(pos, arr, outdir / f"R_{k}_scatter.png", title=f"Residual {k} (scatter)", edges=edge_index)
            if len(pos) >= 3:
                interpolate_to_grid(pos, arr, outdir / f"R_{k}_interp.png", title=f"Residual {k} (interp)")
