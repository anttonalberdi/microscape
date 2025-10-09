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


# ----------------------------- helpers -----------------------------

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


def _index_metabolism_by_spot(meta: Dict[str, Any]) -> Dict[str, Any]:
    """
    Return a dict mapping spot_id-like keys to node dicts with at least {"microbes": {...}}.
    Accepts:
      A) { "S0001": {"microbes": {...}}, ... }
      B) { "spots": { "S0001": {"microbes": {...}}, ... } }
    """
    base = meta.get("spots") if isinstance(meta, dict) and "spots" in meta else meta
    per_spot: Dict[str, Any] = {}
    for k, node in (base or {}).items():
        if not isinstance(node, dict):
            continue
        # normalize
        microbes = (node.get("microbes") or {})
        sid = str(node.get("spot_id") or k)
        per_spot[str(k)] = {"spot_id": sid, "microbes": microbes}
        # also index by explicit spot_id if provided
        per_spot.setdefault(sid, {"spot_id": sid, "microbes": microbes})
    return per_spot


def _spot_join_keys(spot_obj: Dict[str, Any], sid_from_path: str) -> List[str]:
    """
    Build a list of plausible keys to join metabolism JSON by spot.
    Order matters (first hit wins).
    """
    s = spot_obj.get("spot") or spot_obj
    keys: List[str] = []
    # explicit IDs first
    for k in ("id", "name", "spot_id"):
        v = s.get(k)
        if v is not None:
            keys.append(str(v))
    # path stem (S0001)
    keys.append(str(sid_from_path))
    # tolerant variants (strip leading zeros)
    try:
        if sid_from_path.lstrip("S").isdigit():
            n = sid_from_path.lstrip("S")
            keys.append(n)
    except Exception:
        pass
    # unique-ify preserving order
    out: List[str] = []
    seen = set()
    for k in keys:
        if k not in seen:
            out.append(k); seen.add(k)
    return out


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


def _combine_targets(values: List[Optional[float]], op: str, require_all: bool) -> Optional[float]:
    """
    Combine multiple target values with the chosen operator.
    values may include None if some targets are missing.
    """
    if require_all and any(v is None for v in values):
        return None

    vs = [v for v in values if v is not None]
    if not vs:
        return None

    op = (op or "identity").lower()
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


# ----------------------------- CLI -----------------------------

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
    # New: control how we match spot IDs between system/spot files and metabolism JSON
    spot_join: str = typer.Option("auto", help="How to join spots to metabolism JSON: auto|id|name|stem"),
    debug_missing: bool = typer.Option(False, help="Write a CSV of spots/microbes missing targets."),
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
    per_spot = _index_metabolism_by_spot(meta)

    rows: List[Dict[str, Any]] = []
    missing_rows: List[Dict[str, Any]] = []

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
                # Build candidate keys for matching to metabolism JSON
                if spot_join == "auto":
                    keys = _spot_join_keys(spot_obj, sid)
                elif spot_join in ("id", "name", "stem"):
                    if spot_join == "id":
                        keys = [str((spot.get("id") or ""))]
                    elif spot_join == "name":
                        keys = [str((spot.get("name") or ""))]
                    else:
                        keys = [str(sid)]
                else:
                    raise typer.BadParameter("--spot-join must be one of auto|id|name|stem")

                # first existing key wins
                rec_spot = None
                matched_key = None
                for k in keys:
                    if k and k in per_spot:
                        rec_spot = per_spot[k]
                        matched_key = k
                        break

                # environment metadata
                env_meta = {}
                try:
                    env_doc = yaml.safe_load(Path(env_file).read_text())
                    e = env_doc.get("environment", {})
                    env_meta["env_id"] = e.get("id")
                    env_meta["treatment"] = ((e.get("factors") or {}).get("treatment"))
                    env_meta["cage"] = ((e.get("blocking") or {}).get("pen") or (e.get("blocking") or {}).get("cage"))
                    env_meta["weight"] = ((e.get("covariates") or {}).get("weight"))
                except Exception:
                    pass

                if not rec_spot:
                    # No metabolism for this spot; we can still advance progress.
                    if use_progress and progress_ctx is not None:
                        progress_ctx.advance(task_id, 1)
                    continue

                microbes_block = (rec_spot.get("microbes") or {})
                # If metabolism lacks microbes, skip
                if not microbes_block:
                    if use_progress and progress_ctx is not None:
                        progress_ctx.advance(task_id, 1)
                    continue

                # Row assembly
                for mid, mrec in microbes_block.items():
                    vals = [_extract_single_target(mrec, t) for t in targets]
                    combined = _combine_targets(vals, target_op, require_all=combine_require_all)
                    if combined is None:
                        if debug_missing:
                            missing_rows.append({
                                "spot_key": matched_key,
                                "spot_candidates": ";".join(keys),
                                "microbe": mid,
                                "targets": "|".join(targets),
                                "values": json.dumps(vals),
                            })
                        continue

                    feats = _flatten_features_from_spot(
                        spot,
                        microbe_id=mid,
                        include_metabolites=include_metabolites,
                        include_abundance=include_abundance,
                        include_transcripts_sum=include_transcripts_sum,
                    )

                    row = {
                        "spot_id": str(spot.get("name") or spot.get("id") or matched_key or sid),
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

    if debug_missing and missing_rows:
        miss_path = outdir / "missing_targets.csv"
        pd.DataFrame(missing_rows).to_csv(miss_path, index=False)
        typer.secho(f"Wrote debug of missing targets: {miss_path}", fg=typer.colors.YELLOW)

    if not rows:
        # Provide actionable diagnostics
        seen_spots = list(per_spot.keys())[:20]
        typer.secho("No rows assembled. Check your metabolism JSON and target selection.", fg=typer.colors.RED)
        typer.secho(f"First spot keys seen in metabolism JSON: {seen_spots}", fg=typer.colors.BLUE)
        typer.secho("Try: --spot-join stem  or  --spot-join id  depending on how your JSON is keyed.", fg=typer.colors.YELLOW)
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
        "spot_join": spot_join,
    }
    (outdir / "schema.json").write_text(json.dumps(schema, indent=2))

    typer.echo(f"Built modeling table with {len(df)} rows.")
    typer.echo(f"  CSV     : {csv_path}")
    if pq_path:
        typer.echo(f"  Parquet : {pq_path}")
    typer.echo(f"  Schema  : {outdir / 'schema.json'}")