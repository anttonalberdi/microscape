
from __future__ import annotations
import json, subprocess, sys, typer, csv
from pathlib import Path
import typer, numpy as np
from rich.progress import Progress

from ..coupling.loop import run_minimal, compute_summary, save_summary_csv, run_demo_gsmm
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

@app.command("demo2")
def demo2(
    config: str = typer.Option(
        "examples/00_synthetic/community_sbml.yml",
        help="YAML config file for the SBML/dFBA demo",
    ),
    outdir: str = typer.Option(
        "outputs/demo2_gsmm",
        help="Output directory",
    ),
    plot: bool = typer.Option(
        True, help="Save PNG heatmaps if plotting is available"
    ),
):
    """
    Demo 2: SBML‐driven community (GSMM dFBA + diffusion).
    Loads 3 SBML models (FD/LU/BP), runs (p)FBA each step, couples exchanges to RD fields.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    typer.echo("🔧 Initialising Demo 2 (GSMM via SBML)…")
    with Progress() as progress:
        task = progress.add_task("[cyan]Simulating…", total=100)  # coarse bar
        def cb(i, total):
            progress.update(task, completed=min(i, 100))
        fields, summary = run_demo_gsmm(config, outdir, progress_cb=cb)

    # save arrays
    np.savez_compressed(outdir / "fields.npz", **fields)

    # save summary
    with open(outdir / "summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        for k, v in summary.items():
            w.writerow([k, v])

    # optional plots
    if plot:
        try:
            from ..viz.plotting import save_heatmap
            for k in ["glc", "lac", "ac", "but"]:
                if k in fields:
                    save_heatmap(fields[k], outdir / f"{k}.png", title=k.upper())
        except Exception as e:
            typer.echo(f"(plotting skipped: {e})")

    typer.echo("✅ Demo 2 complete.")
    typer.echo(f"   • Arrays:  {outdir/'fields.npz'}")
    typer.echo(f"   • Summary: {outdir/'summary.csv'}")
    if plot:
        typer.echo(f"   • Plots:   {outdir/'glc.png'}, {outdir/'lac.png'}, {outdir/'ac.png'}, {outdir/'but.png'}")