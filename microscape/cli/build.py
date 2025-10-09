from __future__ import annotations
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
import json
import datetime as dt

import typer
import pandas as pd
import yaml
from rich.progress import (
    Progress,
    SpinnerColumn,
    BarColumn,
    TextColumn,
    TimeRemainingColumn,
    MofNCompleteColumn,
)

# Your existing IO helpers
from ..io.system_loader import load_system, iter_spot_files_for_env
from ..io.spot_loader import load_spot


app = typer.Typer(add_completion=False, no_args_is_help=True)


def _flatten_features_from_spot(
    spot: Dict[str, Any],
    microbe_id: str,
    include_metabolites: bool = True,
    include_abundance: bool = True,
    include_transcripts_sum: bool = False,
) -> Dict[str, Any]:
    feats: Dict[str, Any] = {}
    meas = (spot.get("measurements") or {})

    # Microbe abundance (as-is; log1p handled during training)
    if include_abundance:
        mic = (meas.get("microbes") or {}).get("values") or {}
        if isinstance(mic, dict):
            feats["abundance"] = mic.get(microbe_id)

    # Metabolites
    if include_metabolites:
        mets = (meas.get("metabolites") or {}).get("values") or {}
        if isinstance(mets, dict):
            for mid, val in mets.items():
                feats[f"met:{mid}"] = val

    # Simple transcript aggregate (optional)
    if include_transcripts_sum:
        tx = (meas.get("transcripts") or {}).get("values") or {}
        vals = (tx.get(microbe_id) or {})
        if isinstance(vals, dict) and vals:
            try:
                feats["tx_sum"] = float(sum(float(v) for v in vals.values()))
            except Exception:
                feats["tx_sum"] = None

    return feats


def _read_metabolism_json(path: Path) -> Dict[str, Any]:
    try:
        return json.loads(path.read_text())
    except Exception as e:
        raise typer.BadParameter(f"Cannot read metabolism JSON at {path}: {e}")


def _index_metabolism_by_spot(meta: Dict[str, Any]) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    Build two indices:
      - by key:      per_spot[k]
      - by 'spot_id' per_spot_by_sid[ node['spot_id'] or k ]
    Accepts two shapes:
      A) { "S0001": {"microbes": {...}}, ... }
      B) { "spots": { "S0001": {"microbes": {...}}, ... } }
    """
    if "spots" in meta and isinstance(meta["spots"], dict):
        base = {sid: {"microbes": (node or {}).get("microbes", {})} for sid, node in meta["spots"].items()}
    else:
        base = meta

    per_spot: Dict[str, Any] = {}
    for k, node in (base or {}).items():
        if isinstance(node, dict):
            per_spot[str(k)] = node

    per_spot_by_sid: Dict[str, Any] = {}
    for k, node in per_spot.items():
        sid2 = str((node or {}).get("spot_id") or k)
        per_spot_by_sid[sid2] = node

    return per_spot, per_spot_by_sid


def _extract_single_target(mrec: Dict[str, Any], token: str) -> Optional[float]:
    """
    Get a single target value from a microbe record.
      - 'objective' -> rec['objective']
      - 'flux:RID'  -> rec['fluxes'][RID] OR rec['fluxes_all_exchanges'][RID]
    Returns float or None if unavailable.
    """
    if not isinstance(mrec, dict):
        return None

    if token == "objective":
        return mrec.get("objective")

    if token.startswith("flux:"):
        rid = token.split(":", 1)[1]
        fx = (mrec.get("fluxes") or {})
        if rid in fx:
            return fx.get(rid)
        # Fallback to a broader exchange set (if present)
        fx_all = (mrec.get("fluxes_all_exchanges") or {})
        return fx_all.get(rid)

    return None


def _combine_targets(values: List[float], op: str, require_all: bool) -> Optional[float]:
    """
    Combine multiple target values with the chosen operator.
    values may include None if some targets are missing.

    op:
      - identity       -> requires exactly one value
      - sum            -> sum of available values
      - mean           -> mean of available values
      - uptake_sum     -> sum(-v for v<0), i.e., total uptake magnitude (positive)
      - secretion_sum  -> sum(v for v>0), i.e., total secretion (positive)
      - net            -> sum(values), signed

    If require_all=True and any value is None -> return None.
    Otherwise, ignore Nones (if all None -> None).
    """
    if require_all and any(v is None for v in values):
        return None

    vs = [v for v in values if v is not None]
    if not vs:
        return None

    op = op.lower()
    if op == "identity":
        return vs[0] if len(vs) == 1 else None
    if op == "sum":
        return float(sum(vs))
    if op == "mean":
        return float(sum(vs) / len(vs))
    if op == "uptake_sum":
        return float(sum(-v for v in vs if v < 0))
    if op == "secretion_sum":
        return float(sum(v for v in vs if v > 0))
    if op == "net":
        return float(sum(vs))
    # default: identity
    return vs[0] if len(vs) == 1 else None


def _count_spots(env_files: List[str]) -> Optional[int]:
    """Best-effort count of total spots across environment YAMLs (for progress bar total)."""
    total = 0
    try:
        for env_file in env_files:
            env_doc = yaml.safe_load(Path(env_file).read_text())
            env = env_doc.get("environment", {})
            spots = env.get("spots") or []
            total += len(spots)
        return total
    except Exception:
        return None


@app.command("build")
def build_cmd(
    system_yml: Path = typer.Argument(..., help="Path to system.yml"),
    metabolism_json: Path = typer.Option(..., "--metabolism-json", help="Metabolism JSON (un/constrained)."),
    outdir: Path = typer.Option("outputs/model", help="Output folder for the modeling table and schema."),
    targets: List[str] = typer.Option(["objective"], "--target", help="Repeatable. Each is 'objective' or 'flux:EX_*'."),
    target_op: str = typer.Option(
        "identity",
        "--target-op",
        help="How to combine multiple targets: identity|sum|mean|uptake_sum|secretion_sum|net",
    ),
    combine_require_all: bool = typer.Option(
        False, "--combine-require-all/--combine-skip-missing",
        help="If set, drop row unless ALL requested targets are present."
    ),
    include_metabolites: bool = typer.Option(True, help="Include metabolite (met:*) features."),
    include_abundance: bool = typer.Option(True, help="Include microbe abundance feature."),
    include_transcripts_sum: bool = typer.Option(False, help="Include simple transcripts sum per microbe."),
    csv_only: bool = typer.Option(False, help="Write CSV only (skip Parquet)."),
    progress: bool = typer.Option(True, "--progress/--no-progress", help="Show a progress bar while reading spots."),
):
    """
    Build a tidy (spot, microbe) table for modeling by joining spot features
    (metabolites, abundance, transcripts) with per-microbe targets from a metabolism JSON.

    You may pass multiple --target flags and choose how to combine them via --target-op.
    """
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # Load system → environments and paths
    sys_info = load_system(system_yml)
    env_files = sys_info["environment_files"]

    # Load metabolism JSON and index by spot
    meta = _read_metabolism_json(metabolism_json)
    per_spot, per_spot_by_sid = _index_metabolism_by_spot(meta)

    rows: List[Dict[str, Any]] = []

    # Prepare progress bar
    total_spots = _count_spots(env_files)
    use_progress = progress
    progress_ctx = Progress(
        SpinnerColumn(),
        TextColumn("[bold]Building[/bold]"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        MofNCompleteColumn(),
        TimeRemainingColumn(),
        transient=True,
    ) if use_progress else None

    if use_progress and progress_ctx is not None:
        progress_ctx.start()
        task_id = progress_ctx.add_task(
            "Reading spots",
            total=total_spots if isinstance(total_spots, int) and total_spots > 0 else None,
        )

    try:
        for env_file in env_files:
            for sid, spot_path in iter_spot_files_for_env(env_file, sys_info["paths"]):
                spot_obj = load_spot(spot_path) or {}
                spot = spot_obj.get("spot") or spot_obj
                spot_id = str(spot.get("name") or spot.get("id") or sid)

                # environment metadata
                env_meta = {}
                try:
                    env_doc = yaml.safe_load(Path(env_file).read_text())
                    env = env_doc.get("environment", {})
                    env_meta["env_id"] = env.get("id")
                    env_meta["treatment"] = ((env.get("factors") or {}).get("treatment"))
                    env_meta["cage"] = ((env.get("blocking") or {}).get("pen") or (env.get("blocking") or {}).get("cage"))
                    env_meta["weight"] = ((env.get("covariates") or {}).get("weight"))
                except Exception:
                    pass

                # Find metabolism record for this spot
                rec_spot = per_spot.get(spot_id) or per_spot_by_sid.get(spot_id) or {}
                microbes_block = (rec_spot.get("microbes") or {})
                if not microbes_block:
                    if use_progress and progress_ctx is not None:
                        progress_ctx.advance(task_id, 1)
                    continue

                for mid, mrec in microbes_block.items():
                    # Collect requested targets per microbe
                    vals = [_extract_single_target(mrec, t) for t in targets]
                    combined = _combine_targets(vals, target_op, require_all=combine_require_all)
                    if combined is None:
                        continue  # skip if nothing to supervise

                    feats = _flatten_features_from_spot(
                        spot,
                        microbe_id=mid,
                        include_metabolites=include_metabolites,
                        include_abundance=include_abundance,
                        include_transcripts_sum=include_transcripts_sum,
                    )
                    row = {
                        "spot_id": spot_id,
                        "microbe": mid,
                        **env_meta,
                        "target": combined,
                    }
                    # keep raw components (diagnostics)
                    for i, tkn in enumerate(targets):
                        row[f"target_{i+1}"] = vals[i]
                        row[f"target_{i+1}_name"] = tkn

                    row.update(feats)
                    rows.append(row)

                if use_progress and progress_ctx is not None:
                    progress_ctx.advance(task_id, 1)
    finally:
        if use_progress and progress_ctx is not None:
            progress_ctx.stop()

    if not rows:
        seen_spots = sorted(set(per_spot.keys()) | set(per_spot_by_sid.keys()))
        typer.secho(f"No rows assembled. targets={targets}, op={target_op}", fg=typer.colors.RED)
        if seen_spots:
            typer.secho(f"Spots present in metabolism JSON (first 10): {seen_spots[:10]}", fg=typer.colors.YELLOW)
        else:
            typer.secho("No spots were found in the metabolism JSON (unexpected).", fg=typer.colors.YELLOW)
        typer.secho(
            "Tips:\n"
            "  • Try a single '--target objective' first to verify joins.\n"
            "  • For 'flux:EX_*', ensure those EX IDs are recorded in 'fluxes' or 'fluxes_all_exchanges'.\n"
            "  • Use '--combine-skip-missing' (default) to keep rows if at least one requested flux is present,\n"
            "    or '--combine-require-all' to drop rows unless all requested targets exist.",
            fg=typer.colors.BLUE,
        )
        raise typer.Exit(1)

    df = pd.DataFrame(rows)

    # Write CSV
    csv_path = outdir / "table.csv"
    df.to_csv(csv_path, index=False)

    # Parquet (optional)
    pq_path = None
    if not csv_only:
        try:
            pq_path = outdir / "table.parquet"
            df.to_parquet(pq_path, index=False)
        except Exception as e:
            typer.secho(f"Parquet write failed (continuing with CSV): {e}", fg=typer.colors.YELLOW)

    # Schema / provenance
    schema = {
        "created_utc": dt.datetime.utcnow().isoformat(),
        "system": str(system_yml),
        "metabolism_json": str(metabolism_json),
        "targets": targets,
        "target_op": target_op,
        "combine_require_all": combine_require_all,
        "n_rows": int(len(df)),
        "n_features": int(len([c for c in df.columns if c not in {
            "spot_id", "microbe", "env_id", "treatment", "cage", "weight", "target",
            *[f"target_{i+1}" for i in range(len(targets))],
            *[f"target_{i+1}_name" for i in range(len(targets))],
        }])),
        "columns": list(df.columns),
        "types": {c: str(df[c].dtype) for c in df.columns},
    }
    (outdir / "schema.json").write_text(json.dumps(schema, indent=2))

    typer.echo(f"Built modeling table with {len(df)} rows.")
    typer.echo(f"  CSV     : {csv_path}")
    if pq_path:
        typer.echo(f"  Parquet : {pq_path}")
    typer.echo(f"  Schema  : {outdir / 'schema.json'}")
