from __future__ import annotations
import typer, json
from pathlib import Path
import numpy as np
from ..coupling.loop import run_minimal
from ..io import sdp as sdpio

app = typer.Typer(add_completion=False)

@app.command()
def validate_sdp(path: str):
    schema = sdpio.validate_sdp(path)
    typer.echo(json.dumps(schema.model_dump(), indent=2))

@app.command()
def simulate(config: str = typer.Argument(None), out: str = "outputs/run_001.npz"):
    Path(Path(out).parent).mkdir(parents=True, exist_ok=True)
    res = run_minimal()
    np.savez_compressed(out, **res)
    typer.echo(f"Saved results to {out}")

if __name__ == "__main__":
    app()
