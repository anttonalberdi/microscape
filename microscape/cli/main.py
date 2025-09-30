# microscape/cli/main.py
from __future__ import annotations
import json
import subprocess
import sys
from pathlib import Path
import typer
import numpy as np

from ..coupling.loop import run_minimal
from ..io import sdp as sdpio

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
def demo(out: str = "outputs/run_demo.npz"):
    """Run the built-in synthetic demo without any paths."""
    Path(Path(out).parent).mkdir(parents=True, exist_ok=True)
    res = run_minimal()
    np.savez_compressed(out, **res)
    typer.echo(f"Demo finished. Saved results to {out}")

@app.command()
def update(
    ref: str = typer.Option("main", help="Git ref: branch, tag, or commit."),
    with_deps: bool = typer.Option(False, help="Also reinstall dependencies."),
    yes: bool = typer.Option(False, "--yes", "-y", help="Skip confirmation."),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Less pip output."),
    no_build_iso: bool = typer.Option(True, help="Use --no-build-isolation for speed."),
):
    """
    Update MicroScape from GitHub into the current environment.
    By default, only the package is reinstalled (no dependency reinstall).
    """
    url = f"{REPO_URL}@{ref}"
    if not yes:
        typer.confirm(
            f"This will install from {url} into the current environment ({sys.executable}). Continue?",
            abort=True,
        )

    cmd = [sys.executable, "-m", "pip", "install", "--upgrade", "--force-reinstall", url]
    if not with_deps:
        cmd.insert(5, "--no-deps")  # after 'install'
    if no_build_iso:
        cmd.insert(5, "--no-build-isolation")
    if quiet:
        cmd.insert(5, "-q")

    typer.echo("Running: " + " ".join(cmd))
    try:
        subprocess.check_call(cmd)
        typer.secho("MicroScape updated successfully.", fg=typer.colors.GREEN)
    except subprocess.CalledProcessError as e:
        typer.secho(f"Update failed (exit {e.returncode}).", fg=typer.colors.RED)
        raise typer.Exit(code=e.returncode)