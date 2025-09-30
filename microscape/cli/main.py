from __future__ import annotations
import typer, numpy as np
from pathlib import Path
from ..coupling.loop import run_minimal
from ..utils import get_demo_dir

app = typer.Typer(add_completion=False)

@app.command()
def demo(out: str = "outputs/run_demo.npz"):
    """
    Run the bundled synthetic demo without specifying any paths.
    """
    demo_path = get_demo_dir()
    typer.echo(f"Using demo at: {demo_path}")
    Path(Path(out).parent).mkdir(parents=True, exist_ok=True)
    res = run_minimal()
    np.savez_compressed(out, **res)
    typer.echo(f"Saved results to {out}")
