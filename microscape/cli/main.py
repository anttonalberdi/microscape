
from __future__ import annotations
import json, subprocess, sys, typer
from pathlib import Path
import typer, numpy as np
from rich.progress import Progress

from ..coupling.loop import run_minimal, compute_summary, save_summary_csv
from ..io import sdp as sdpio
from ..viz.plotting import save_heatmap, save_profile

app = typer.Typer(add_completion=False)
REPO_URL = "git+https://github.com/anttonalberdi/microscape.git"

@app.command()
def validate_sdp(path: str):
    """Validate a Spatial Data Package (SDP)."""
    schema = sdpio.validate_sdp(path)
    typer.echo(json.dumps(schema.model_dump(), indent=2))

@app.command()
def simulate(config: str = typer.Argument(None), out: str = "outputs/run_001.npz"):
    """Run a minimal synthetic simulation demo (placeholder for config-driven runs)."""
    Path(Path(out).parent).mkdir(parents=True, exist_ok=True)
    res = run_minimal()
    np.savez_compressed(out, **res)
    typer.echo(f"Saved results to {out}")

@app.command()
def demo(outdir: str = "outputs/demo_001", plot: bool = typer.Option(True, help="Save PNG plots")):
    """
    Run the built-in synthetic demo and produce user-friendly outputs:
    - fields.npz            (butyrate, fibre, mucosa_mask)
    - summary.csv           (means, SCFA at mucosa)
    - butyrate.png          (heatmap)
    - radial_profile.png    (butyrate vs distance from mucosa)
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    typer.echo("🔧 Initialising demo (fibre → butyrate on a small tissue slice)…")

    # Progress bar for the simulation loop
    steps = 200
    with Progress() as progress:
        task = progress.add_task("[cyan]Simulating diffusion & production…", total=steps)
        def cb(i, total):
            progress.update(task, completed=i)
        res = run_minimal(sim_steps=steps, dt=5.0, voxel_um=10.0, progress_cb=cb)

    typer.echo("🧮 Computing summaries (SCFA-at-mucosa, radial profile)…")
    summary = compute_summary(res["butyrate"], res["mucosa_mask"])

    typer.echo("💾 Writing outputs (arrays, CSV, plots)…")
    fields_npz = outdir / "fields.npz"
    np.savez_compressed(fields_npz, **res)
    save_summary_csv(summary, outdir / "summary.csv")

    if plot:
        save_heatmap(res["butyrate"], outdir / "butyrate.png", title="Butyrate (a.u.)")
        save_heatmap(res["fibre"], outdir / "fibre.png", title="Fibre", cmap="magma")
        save_profile(summary["profile_dist_px"], summary["profile_mean"], outdir / "radial_profile.png",
                     title="Butyrate vs distance from mucosa")

    typer.echo("✅ Demo complete.")
    typer.echo(f"   • Arrays:   {fields_npz}")
    typer.echo(f"   • Summary:  {outdir/'summary.csv'}")
    if plot:
        typer.echo(f"   • Plots:    {outdir/'butyrate.png'}, {outdir/'radial_profile.png'}")

@app.command()
def update(
    ref: str = typer.Option("main", help="Git ref: branch, tag, or commit."),
    yes: bool = typer.Option(False, "--yes", "-y", help="Skip confirmation."),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Less pip output."),
    with_deps: bool = typer.Option(False, help="Reinstall dependencies too."),
):
    url = f"{REPO_URL}@{ref}"
    if not yes:
        typer.confirm(f"This will install from {url} into ({sys.executable}). Continue?", abort=True)

    base = [sys.executable, "-m", "pip", "install", "--upgrade", "--force-reinstall"]
    if quiet:
        base.append("-q")

    # 1) Fast path: no deps, no build isolation
    cmd = base + ["--no-deps", "--no-build-isolation", url]
    typer.echo("Running: " + " ".join(cmd))
    try:
        subprocess.check_call(cmd)
        typer.secho("MicroScape updated (fast path).", fg=typer.colors.GREEN)
        return
    except subprocess.CalledProcessError:
        typer.secho("Fast update failed; retrying with build isolation…", fg=typer.colors.YELLOW)

    # 2) Retry with build isolation (installs hatchling if needed); include deps if requested
    cmd = base + ([url])
    typer.echo("Running: " + " ".join(cmd))
    try:
        subprocess.check_call(cmd)
        typer.secho("MicroScape updated (build isolation).", fg=typer.colors.GREEN)
    except subprocess.CalledProcessError as e:
        typer.secho(
            f"Update failed (exit {e.returncode}). "
            "Install 'hatchling' in your env or run with '--with-deps'.",
            fg=typer.colors.RED,
        )
        raise typer.Exit(code=e.returncode)