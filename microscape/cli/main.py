
from __future__ import annotations
import json, subprocess, sys, typer, csv, shutil
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
        help="YAML config for the SBML/dFBA demo",
    ),
    outdir: str = typer.Option("outputs/demo2_gsmm", help="Output directory"),
    plot: bool = typer.Option(True, help="Save PNG heatmaps"),
):
    # --- lazy import so the rest of the CLI works even if COBRA or the file is missing ---
    try:
        # IMPORTANT: import here, not at top-level
        from ..coupling.loop_gsmm import run_demo_gsmm
    except Exception as e:
        typer.secho(
            "GSMM demo not available. Ensure file 'microscape/coupling/loop_gsmm.py' "
            "is packaged and COBRApy is installed (pip install cobra optlang).",
            fg=typer.colors.RED,
        )
        raise typer.Exit(code=1)

    from rich.progress import Progress
    from pathlib import Path
    import numpy as np, csv

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    typer.echo("🔧 Initialising Demo 2 (GSMM via SBML)…")
    with Progress() as progress:
        task = progress.add_task("[cyan]Simulating…", total=100)
        def cb(i, total): progress.update(task, completed=min(i, 100))
        fields, summary = run_demo_gsmm(config, outdir, progress_cb=cb)

    np.savez_compressed(outdir / "fields.npz", **fields)
    with open(outdir / "summary.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["metric","value"]); [w.writerow([k,v]) for k,v in summary.items()]

    if plot:
        try:
            from ..viz.plotting import save_heatmap
            for k in ["glc","lac","ac","but"]:
                if k in fields: save_heatmap(fields[k], outdir / f"{k}.png", title=k.upper())
        except Exception as e:
            typer.echo(f"(plotting skipped: {e})")

    typer.echo("✅ Demo 2 complete.")
