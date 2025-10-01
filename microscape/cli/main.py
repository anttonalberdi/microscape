
from __future__ import annotations
from pathlib import Path
import json
import typer
from rich.progress import Progress
import numpy as np

from ..core.registry import get as get_engine
from ..kinetics import toy_chain  # ensure "toy" is registered
from ..runner.simulate_graph import simulate_graph
from ..io.graph_config import load_graph_yaml
from ..viz.graph import scatter_field, interpolate_to_grid

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

@app.command("simulate")
def simulate_cmd(
    config: Path = typer.Argument(..., help="YAML config (graph or grid)"),
    engine: str = typer.Option("toy", help="Engine: toy (others can be added later)"),
    outdir: Path = typer.Option("outputs/run", help="Output directory"),
    plot: bool = typer.Option(True, help="Save plots"),
):
    """
    General simulation entrypoint. Dispatches on space.type in the YAML.
    """
    # Load once (also lets us resolve relative paths in future)
    cfg = load_graph_yaml(config)  # works for graph; harmless for reading 'space' key for grid later
    space_type = (cfg.get("space") or {}).get("type", "graph").lower()

    step = get_engine(engine)
    outdir.mkdir(parents=True, exist_ok=True)

    with Progress() as progress:
        task = progress.add_task(f"[cyan]Simulating ({space_type})…", total=100)
        cb = lambda i, total: progress.update(task, completed=min(i, 100))

        if space_type == "graph":
            fields, summary = simulate_graph(str(config), step, progress=cb)
        else:
            raise typer.BadParameter(
                f"Unsupported space.type '{space_type}'. "
                "Currently supported: 'graph'."
            )

    # Save outputs
    np.savez_compressed(outdir / "fields.npz", **fields)
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2))

    # Optional plotting (graph)
    if plot and space_type == "graph":
        nodes = cfg["space"]["nodes"]; edges = cfg["space"]["edges"]
        pos = np.array([nd.get("pos_um", [0, 0])[:2] for nd in nodes], float)
        id_to_idx = {nd["id"]: i for i, nd in enumerate(nodes)}
        edge_index = np.array([[id_to_idx[e["i"]], id_to_idx[e["j"]]] for e in edges], int)
        for k, arr in fields.items():
            scatter_field(pos, arr, outdir / f"{k}_scatter.png", title=f"{k} (scatter)", edges=edge_index)
            interpolate_to_grid(pos, arr, outdir / f"{k}_interp.png", title=f"{k} (interp)")

    typer.secho("✅ Done.", fg=typer.colors.GREEN)