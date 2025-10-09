
from __future__ import annotations
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple

import json
import numpy as np
import pandas as pd
import typer

from rich.progress import (
    Progress,
    SpinnerColumn,
    BarColumn,
    TextColumn,
    TimeRemainingColumn,
    MofNCompleteColumn,
)

# Reuse predict-time helpers (local copies to keep module standalone)
def _standardize_with_stats(X: pd.DataFrame, means: Dict[str, float], sds: Dict[str, float]) -> pd.DataFrame:
    Z = X.copy()
    for c in Z.columns:
        mu = float(means.get(c, Z[c].mean()))
        sd = float(sds.get(c, Z[c].std(ddof=0) or 1.0))
        if sd == 0:
            sd = 1.0
        Z[c] = (Z[c].astype(float) - mu) / sd
    return Z

def _rff_apply(coords: np.ndarray, params: Dict[str, Any]) -> pd.DataFrame:
    D = int(params["D"])
    omega = np.array(params["omega"], dtype=float)  # (D, d)
    b = np.array(params["b"], dtype=float)          # (D,)
    proj = coords @ omega.T + b
    phi = np.sqrt(2.0 / D) * np.cos(proj)
    cols = [f"rff:{i}" for i in range(D)]
    return pd.DataFrame(phi, columns=cols)

def _posterior_array(idata, name: str):
    try:
        return idata.posterior[name].to_numpy()
    except Exception:
        return None

def _encode_with_map(series: pd.Series, mapping: Dict[str, int]):
    if not mapping: return None
    return series.fillna("__NA__").astype(str).map(mapping).values.astype(int)

def _hdi(x: np.ndarray, prob: float = 0.95, axis: Optional[int] = None):
    lo = (1.0 - prob) / 2.0
    hi = 1.0 - lo
    q = np.quantile(x, [lo, hi], axis=axis)
    return q[0], q[1]


app = typer.Typer(add_completion=False, no_args_is_help=True)


def _load_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path)


def _apply_scope_mask(df: pd.DataFrame, scope: str) -> np.ndarray:
    """scope: 'all' | 'spot:S0001' | 'env:E001'"""
    if scope == "all" or not scope:
        return np.ones(len(df), dtype=bool)
    if scope.startswith("spot:"):
        sid = scope.split(":", 1)[1]
        return (df.get("spot_id") == sid).values
    if scope.startswith("env:"):
        eid = scope.split(":", 1)[1]
        return (df.get("env_id") == eid).values
    return np.ones(len(df), dtype=bool)


@app.command("whatif")
def whatif_cmd(
    table: Path = typer.Argument(..., help="Existing modeling table (CSV or Parquet)."),
    outdir: Path = typer.Option("outputs/whatif", help="Where to save the modified table and (optionally) predictions."),
    # Batch operations
    drop_microbe: List[str] = typer.Option(None, "--drop-microbe", help="Remove a microbe across scope (rows where microbe==ID)."),
    zero_metabolite: List[str] = typer.Option(None, "--zero-metabolite", help="Set a metabolite feature 'met:ID' to 0 across scope."),
    set_abundance: List[str] = typer.Option(None, "--set-abundance", help="microbe=value; set abundance for that microbe within scope."),
    scale_abundance: List[str] = typer.Option(None, "--scale-abundance", help="microbe=factor; multiply abundance within scope."),
    set_metabolite: List[str] = typer.Option(None, "--set-metabolite", help="Cxxxx=value; set 'met:Cxxxx' within scope."),
    scale_metabolite: List[str] = typer.Option(None, "--scale-metabolite", help="Cxxxx=factor; scale 'met:Cxxxx' within scope."),
    scope: str = typer.Option("all", help="Scope: all | spot:S0001 | env:E001"),
    # Optional prediction in one go
    posterior_nc: Optional[Path] = typer.Option(None, help="If provided, run prediction after edits."),
    model_card: Optional[Path] = typer.Option(None, help="Required with --posterior-nc to reconstruct features/RFF."),
    hdi_prob: float = typer.Option(0.95, help="HDI for predictions."),
    include_draws: bool = typer.Option(False, help="Write per-draw predictions too (large)."),
    progress: bool = typer.Option(True, "--progress/--no-progress", help="Show progress bar."),
):
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    df = _load_table(table).copy()
    mask = _apply_scope_mask(df, scope)

    # --- Apply edits ---
    if drop_microbe:
        before = len(df)
        df = df[~(df["microbe"].isin(drop_microbe) & mask)].copy()
        typer.secho(f"Dropped {before - len(df)} row(s) matching microbes in scope.", fg=typer.colors.YELLOW)

    def _parse_pairs(kvs: List[str], kind: str) -> Dict[str, float]:
        parsed = {}
        for kv in kvs or []:
            if "=" not in kv:
                raise typer.BadParameter(f"--{kind} expects entries like name=value, got '{kv}'")
            k, v = kv.split("=", 1)
            parsed[k.strip()] = float(v)
        return parsed

    # abundance is per (spot,microbe)
    set_ab = _parse_pairs(set_abundance, "set-abundance")
    sc_ab = _parse_pairs(scale_abundance, "scale-abundance")

    if set_ab:
        for m, val in set_ab.items():
            sel = (df["microbe"] == m) & mask
            df.loc[sel, "abundance"] = float(val)
    if sc_ab:
        for m, fac in sc_ab.items():
            sel = (df["microbe"] == m) & mask
            df.loc[sel, "abundance"] = df.loc[sel, "abundance"].astype(float) * float(fac)

    # metabolites are spot-level columns 'met:*'
    set_met = _parse_pairs(set_metabolite, "set-metabolite")
    sc_met = _parse_pairs(scale_metabolite, "scale-metabolite")

    for cid, val in set_met.items():
        col = f"met:{cid}" if not cid.startswith("met:") else cid
        if col in df.columns:
            df.loc[mask, col] = float(val)
    for cid, fac in sc_met.items():
        col = f"met:{cid}" if not cid.startswith("met:") else cid
        if col in df.columns:
            df.loc[mask, col] = df.loc[mask, col].astype(float) * float(fac)

    for cid in zero_metabolite or []:
        col = f"met:{cid}" if not cid.startswith("met:") else cid
        if col in df.columns:
            df.loc[mask, col] = 0.0

    # Save modified table (CSV + Parquet if possible)
    csv_path = outdir / "table_modified.csv"
    df.to_csv(csv_path, index=False)
    try:
        pq_path = outdir / "table_modified.parquet"
        df.to_parquet(pq_path, index=False)
    except Exception:
        pq_path = None

    typer.echo(f"Modified table saved: {csv_path}")
    if pq_path:
        typer.echo(f"Parquet: {pq_path}")

    # --- Optional prediction step (uses spatial RFF if model_card contains it) ---
    if posterior_nc and model_card:
        card = json.loads(Path(model_card).read_text())
        features: List[str] = card.get("features") or []
        feature_means = card.get("feature_means") or {}
        feature_sds = card.get("feature_sds") or {}
        spatial = card.get("spatial") or {}
        group_maps = card.get("group_maps") or {}
        cage_map = group_maps.get("cage") or {}
        env_map = group_maps.get("env") or {}
        treat_map = group_maps.get("treatment") or {}

        # Rebuild X like in predict
        base_cols = [c for c in features if not c.startswith("rff:")]
        rff_cols = [c for c in features if c.startswith("rff:")]

        # Ensure base columns present
        missing = [c for c in base_cols if c not in df.columns]
        if missing:
            raise typer.BadParameter(f"Modified table missing required base features: {missing}")

        X_base = df[base_cols].copy()
        if "abundance" in X_base.columns:
            if (X_base["abundance"] < -1).any():
                raise typer.BadParameter("'abundance' has values < -1; log1p undefined.")
            X_base["abundance"] = np.log1p(X_base["abundance"].astype(float))
        X_base = _standardize_with_stats(X_base, feature_means, feature_sds)

        if rff_cols:
            rff = spatial or {}
            if not (rff.get("rff_gp") and rff.get("rff_params")):
                raise typer.BadParameter("Model used rff:* features but rff_params missing.")
            rff_params = rff["rff_params"]
            coord_cols = rff_params.get("coord_cols") or [c for c in ["x_um","y_um","z_um"] if c in df.columns]
            if len(coord_cols) < 2:
                raise typer.BadParameter("RFF model expects at least x_um and y_um in the table.")
            coords = df[coord_cols].to_numpy(dtype=float)
            Phi = _rff_apply(coords, rff_params)
            Phi = Phi[rff_cols]
            X = pd.concat([X_base, Phi], axis=1)
            X = X[features]
        else:
            X = X_base

        # Encodings
        cage_idx = _encode_with_map(df["cage"], cage_map) if "cage" in df.columns and cage_map else None
        env_idx  = _encode_with_map(df["env_id"], env_map) if "env_id" in df.columns and env_map else None
        tr_idx   = _encode_with_map(df["treatment"], treat_map) if "treatment" in df.columns and treat_map else None

        try:
            import arviz as az
        except Exception as e:
            raise typer.BadParameter(f"arviz is required for prediction: {e}")
        idata = az.from_netcdf(posterior_nc)

        alpha = _posterior_array(idata, "alpha")
        beta  = _posterior_array(idata, "beta")
        a_cage = _posterior_array(idata, "a_cage")
        a_env  = _posterior_array(idata, "a_env")
        gamma_t = _posterior_array(idata, "gamma_t")

        if alpha is None or beta is None:
            raise typer.BadParameter("Posterior missing 'alpha' or 'beta'.")

        n_obs = df.shape[0]
        n_feat = beta.shape[-1]
        if X.shape[1] != n_feat:
            raise typer.BadParameter(f"Feature mismatch: X has {X.shape[1]}, beta has {n_feat}.")

        use_progress = progress
        prog = Progress(
            SpinnerColumn(),
            TextColumn("[bold]What-if predict[/bold]"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            MofNCompleteColumn(),
            TimeRemainingColumn(),
            transient=True,
        ) if use_progress else None

        if use_progress and prog is not None:
            prog.start()
            task = prog.add_task("Computing posterior predictions", total=n_obs)

        chunk = max(1, n_obs // 20)
        mu_chunks = []
        for start in range(0, n_obs, chunk):
            end = min(n_obs, start + chunk)
            X_chunk = X.values[start:end, :]
            base = alpha[..., None] + np.einsum("cdk,mk->cdm", beta, X_chunk)
            if a_cage is not None and cage_idx is not None:
                idx = cage_idx[start:end]
                base = base + a_cage[:, :, idx]
            if a_env is not None and env_idx is not None:
                idx = env_idx[start:end]
                base = base + a_env[:, :, idx]
            if gamma_t is not None and tr_idx is not None:
                idx = tr_idx[start:end]
                base = base + gamma_t[:, :, idx]
            mu_chunks.append(base)
            if use_progress and prog is not None:
                prog.advance(task, end - start)

        if use_progress and prog is not None:
            prog.stop()

        mu = np.concatenate(mu_chunks, axis=-1)  # (C,D,N)
        flat = mu.reshape(-1, n_obs)
        mean = flat.mean(axis=0)
        lo, hi = _hdi(flat, prob=hdi_prob, axis=0)

        meta_cols = [c for c in ["spot_id", "microbe", "env_id", "cage", "treatment"] if c in df.columns]
        out_df = df[meta_cols].copy() if meta_cols else pd.DataFrame(index=df.index)
        out_df["y_hat_mean"] = mean
        out_df[f"y_hat_hdi{int(hdi_prob*100)}_low"] = lo
        out_df[f"y_hat_hdi{int(hdi_prob*100)}_high"] = hi

        pred_path = outdir / "whatif_predictions.csv"
        out_df.to_csv(pred_path, index=False)
        typer.echo(f"Predictions written: {pred_path}")

        if include_draws:
            npz_path = outdir / "whatif_predictions_draws.npz"
            np.savez_compressed(npz_path, draws=flat.T)
            typer.echo(f"Per-draw matrix: {npz_path}")

    else:
        typer.echo("No posterior/model_card given; skipped prediction.")


if __name__ == "__main__":
    app()
