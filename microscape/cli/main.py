from __future__ import annotations
from pathlib import Path
import json, typer, shutil, sys, subprocess
import numpy as np
from rich.progress import Progress
from ..runner.snapshot import run_snapshot
from ..io.graph_config import load_graph_yaml
from ..viz.graph import scatter_field, interpolate_to_grid
from ..validation.system import validate_system
from ..profile.ecology import profile_system

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
    ecology_config: Path = typer.Option(None, "--ecology-config", "-e",
        help="(Optional) override for ecology rules YAML. If omitted, uses system.config.ecology."),
    outdir: Path = typer.Option("work/profile", help="Output directory for enriched spot YAMLs and summaries."),
    overwrite: bool = typer.Option(False, help="Overwrite enriched spot files if they already exist."),
    json_out: bool = typer.Option(True, help="Also write a JSON summary."),
):
    """
    Run ecology profiling:
      - Reads system.yml -> environments -> spots
      - Resolves ecology rules from system.config.ecology (or --ecology-config override)
      - Infers per-microbe trait states/scores from transcripts
      - Writes enriched spots under OUTDIR and a summary CSV/JSON
    """
    system_yml = system_yml.resolve()

    # Resolve ecology rules path
    if ecology_config is None:
        # Read system.yml to find system.config.ecology
        try:
            import yaml
            sysd = yaml.safe_load(system_yml.read_text())
            ec_rel = (((sysd or {}).get("system") or {}).get("config") or {}).get("ecology")
        except Exception as e:
            typer.secho(f"Failed to read {system_yml}: {e}", fg=typer.colors.RED)
            raise typer.Exit(code=2)

        if ec_rel is None:
            typer.secho(
                "No ecology rules provided: system.config.ecology missing and no --ecology-config given.",
                fg=typer.colors.RED
            )
            raise typer.Exit(code=2)

        ecology_config = (system_yml.parent / ec_rel).resolve()

    if not ecology_config.exists():
        typer.secho(f"Ecology config not found: {ecology_config}", fg=typer.colors.RED)
        raise typer.Exit(code=2)

    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    typer.secho(
        f"🧭 Profiling ecology\n  system : {system_yml}\n  rules  : {ecology_config}\n  out    : {outdir}",
        fg=typer.colors.CYAN
    )

    # Optional: quick validation to catch broken references before profiling
    try:
        _, errors, _ = validate_system(system_yml)
        if errors:
            typer.secho("Validation errors found; profiling may fail:", fg=typer.colors.YELLOW)
            for e in errors:
                typer.echo(f"  - {e}")
    except Exception:
        # validation is optional; don’t block profiling
        pass

    with Progress() as progress:
        task = progress.add_task("[cyan]Profiling…", total=100)
        result = profile_system(
            system_yml,
            ecology_config,
            outdir,
            overwrite=overwrite,
            progress_cb=lambda p: progress.update(task, completed=min(100, p)),
        )

    # Write summaries
    (outdir / "profile_summary.csv").write_text(result["summary_csv"])
    if json_out:
        (outdir / "profile_summary.json").write_text(json.dumps(result["summary_json"], indent=2))

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
