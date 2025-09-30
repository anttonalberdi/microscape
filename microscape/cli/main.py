from __future__ import annotations
import json, subprocess, sys
from pathlib import Path
import typer, numpy as np

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

    res = run_minimal()
    fields_npz = outdir / "fields.npz"
    np.savez_compressed(fields_npz, **res)

    summary = compute_summary(res["butyrate"], res["mucosa_mask"])
    save_summary_csv(summary, outdir / "summary.csv")

    if plot:
        save_heatmap(res["butyrate"], outdir / "butyrate.png", title="Butyrate (a.u.)")
        save_heatmap(res["fibre"], outdir / "fibre.png", title="Fibre", cmap="magma")
        save_profile(summary["profile_dist_px"], summary["profile_mean"], outdir / "radial_profile.png",
                     title="Butyrate vs distance from mucosa")

    typer.echo(f"Demo complete.\n- Arrays: {fields_npz}\n- Summary: {outdir/'summary.csv'}\n- Plots: {outdir/'butyrate.png'}, {outdir/'radial_profile.png'}")

@app.command()
def update(
    ref: str = typer.Option("main", help="Git ref: branch, tag, or commit."),
    yes: bool = typer.Option(False, "--yes", "-y", help="Skip confirmation."),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Less pip output."),
    with_deps: bool = typer.Option(False, help="Reinstall dependencies too."),
):
    """
    Update MicroScape from GitHub. Fast path (no deps) first; auto-retry with build isolation if needed.
    """
    url = f"{REPO_URL}@{ref}"
    if not yes:
        typer.confirm(f"This will install from {url} into ({sys.executable}). Continue?", abort=True)

    base = [sys.executable, "-m", "pip", "install", "--upgrade", "--force-reinstall"]
    if quiet:
        base.append("-q")

    # 1) Fast path: no deps, no build isolation
    cmd = base + ["--no-deps", "--no-build-isolation", url]
    try:
        _run(cmd)
        typer.secho("MicroScape updated (fast path).", fg=typer.colors.GREEN)
        return
    except subprocess.CalledProcessError as e:
        errmsg = str(e)
        typer.secho("Fast update failed; retrying with build isolation…", fg=typer.colors.YELLOW)

    # 2) Retry: allow build isolation (will fetch hatchling); include deps if requested
    cmd = base + ([url] if not with_deps else [url])
    try:
        _run(cmd)
        typer.secho("MicroScape updated (build isolation).", fg=typer.colors.GREEN)
    except subprocess.CalledProcessError as e:
        typer.secho(f"Update failed (exit {e.returncode}). Try installing 'hatchling' in your env.", fg=typer.colors.RED)
        raise typer.Exit(code=e.returncode)