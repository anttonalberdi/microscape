# microscape/cli/constrain.py
from __future__ import annotations
from pathlib import Path
import typer
from rich.progress import Progress

from ..io.system_loader import load_system
from ..io.metabolism_rules import load_rules
from ..runner.constraints import constrain_one

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command("constrain")
def constrain_cmd(
    system_yml: Path = typer.Argument(..., help="Path to system.yml"),
    mode: str = typer.Option("environmental", help="environmental | transcriptional | combined"),
    outdir: Path = typer.Option("outputs/constraints", help="Output directory"),
    write_models: bool = typer.Option(False, help="Write constrained SBMLs"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Extra logs"),
):
    """
    Apply constraints per spot×microbe.

    Modes:
      • environmental  – set EX bounds from metabolite mM via metabolism.yml
      • transcriptional – scale internal reaction bounds by gene expression (GPR-aware)
      • combined       – apply both environmental and transcriptional
    """
    sys_info = load_system(system_yml)
    metab_cfg = (sys_info["system"].get("config") or {}).get("metabolism")
    if not metab_cfg:
        typer.secho("No metabolism config set in system.config.metabolism", fg=typer.colors.RED)
        raise typer.Exit(1)
    metab_path = (sys_info["root"] / sys_info["paths"].get("config_dir", "config") / metab_cfg).resolve()
    rules = load_rules(metab_path)

    outdir.mkdir(parents=True, exist_ok=True)
    out_csv = outdir / f"constraints__{mode}.tsv"

    with Progress() as prog:
        task = prog.add_task(f"[cyan]Constraining models ({mode})…", total=100)
        _ = constrain_one(system_yml, rules, out_csv, mode=mode, write_models=write_models, verbose=verbose)
        prog.update(task, completed=100)

    typer.secho("✅ Constraint profiling complete.", fg=typer.colors.GREEN)
    typer.echo(f"  TSV  : {out_csv}")
    if write_models:
        typer.echo(f"  SBML : {outdir / 'models_constrained'}")
