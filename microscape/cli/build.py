from __future__ import annotations
from pathlib import Path
from typing import Dict, Any, List, Optional
import json, typer, datetime as dt
import pandas as pd

from ..io.system_loader import load_system, iter_spot_files_for_env
from ..io.spot_loader import load_spot

app = typer.Typer(add_completion=False, no_args_is_help=True)

def _flatten_features_from_spot(spot: Dict[str, Any], microbe_id: str, include_metabolites: bool=True,
                                include_abundance: bool=True, include_transcripts_sum: bool=False) -> Dict[str, Any]:
    feats: Dict[str, Any] = {}
    meas = (spot.get("measurements") or {})
    if include_abundance:
        mic = (meas.get("microbes") or {}).get("values") or {}
        if isinstance(mic, dict):
            feats["abundance"] = mic.get(microbe_id)
    if include_metabolites:
        mets = (meas.get("metabolites") or {}).get("values") or {}
        if isinstance(mets, dict):
            for mid, val in mets.items():
                feats[f"met:{mid}"] = val
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

def _extract_target(rec: Dict[str, Any], target: str) -> Optional[float]:
    if target == "objective":
        return rec.get("objective", None)
    if target.startswith("flux:"):
        rid = target.split(":",1)[1]
        return (rec.get("fluxes") or {}).get(rid, None)
    return None

@app.command("build")
def build_cmd(
    system_yml: Path = typer.Argument(..., help="Path to system.yml"),
    metabolism_json: Path = typer.Option(..., "--metabolism-json", help="Named metabolism JSON (un/constrained)."),
    outdir: Path = typer.Option("outputs/model", help="Output folder for modeling table."),
    target: str = typer.Option("objective", help="Target variable: 'objective' or 'flux:EX_*'."),
    include_metabolites: bool = typer.Option(True, help="Include metabolites (met:*) features."),
    include_abundance: bool = typer.Option(True, help="Include microbe abundance feature."),
    include_transcripts_sum: bool = typer.Option(False, help="Include simple transcripts sum per microbe."),
    csv_only: bool = typer.Option(False, help="Write CSV only (skip Parquet)."),
):
    """Build a tidy table for modeling by joining spot features with per-microbe targets."""
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    sys_info = load_system(system_yml)
    env_files = sys_info["environment_files"]

    meta = _read_metabolism_json(metabolism_json)
    per_spot = meta
    if "spots" in meta:
        per_spot = {sid: {"microbes": m.get("microbes", {})} for sid, m in meta["spots"].items()}

    rows: List[Dict[str, Any]] = []
    for env_file in env_files:
        for sid, spot_path in iter_spot_files_for_env(env_file, sys_info["paths"]):
            spot_obj = load_spot(spot_path) or {}
            spot = spot_obj.get("spot") or spot_obj
            spot_id = spot.get("name") or spot.get("id") or sid
            # env meta
            env_meta = {}
            try:
                import yaml
                env_doc = yaml.safe_load(Path(env_file).read_text())
                env = env_doc.get("environment", {})
                env_meta["env_id"] = env.get("id")
                env_meta["treatment"] = ((env.get("factors") or {}).get("treatment"))
                env_meta["cage"] = ((env.get("blocking") or {}).get("pen") or (env.get("blocking") or {}).get("cage"))
                env_meta["weight"] = ((env.get("covariates") or {}).get("weight"))
            except Exception:
                pass

            microbes_block = (per_spot.get(spot_id) or {}).get("microbes") or {}
            if not microbes_block:
                continue
            for mid, mrec in microbes_block.items():
                target_val = _extract_target(mrec, target)
                if target_val is None:
                    continue
                feats = _flatten_features_from_spot(spot, microbe_id=mid,
                                                    include_metabolites=include_metabolites,
                                                    include_abundance=include_abundance,
                                                    include_transcripts_sum=include_transcripts_sum)
                row = {
                    "spot_id": spot_id,
                    "microbe": mid,
                    **env_meta,
                    "target": target_val,
                    **feats,
                }
                rows.append(row)

    if not rows:
        typer.secho("No rows assembled. Check your metabolism JSON and target selection.", fg=typer.colors.RED)
        raise typer.Exit(1)

    df = pd.DataFrame(rows)

    csv_path = outdir / "table.csv"
    df.to_csv(csv_path, index=False)

    pq_path = None
    if not csv_only:
        try:
            pq_path = outdir / "table.parquet"
            df.to_parquet(pq_path, index=False)
        except Exception as e:
            typer.secho(f"Parquet write failed (will continue with CSV): {e}", fg=typer.colors.YELLOW)

    schema = {
        "created_utc": dt.datetime.utcnow().isoformat(),
        "system": str(system_yml),
        "metabolism_json": str(metabolism_json),
        "target": target,
        "n_rows": len(df),
        "n_features": len([c for c in df.columns if c not in {"spot_id","microbe","env_id","treatment","cage","weight","target"}]),
        "columns": list(df.columns),
        "types": {c: str(df[c].dtype) for c in df.columns},
    }
    (outdir / "schema.json").write_text(json.dumps(schema, indent=2))

    typer.echo(f"Built modeling table with {len(df)} rows.")
    typer.echo(f"  CSV   : {csv_path}")
    if pq_path:
        typer.echo(f"  Parquet: {pq_path}")
    typer.echo(f"  Schema: {outdir / 'schema.json'}")
