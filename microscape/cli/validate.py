from __future__ import annotations
from pathlib import Path
import json, typer
from ..validate.system import validate_system
from rich.progress import Progress

app = typer.Typer(add_completion=False)

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
