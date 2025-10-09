
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

app = typer.Typer(add_completion=False, no_args_is_help=True)


# ---------------- utils ----------------

def _hdi(x: np.ndarray, prob: float = 0.95, axis: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
    lo = (1.0 - prob) / 2.0
    hi = 1.0 - lo
    q = np.quantile(x, [lo, hi], axis=axis)
    return q[0], q[1]


def _standardize_with_stats(X: pd.DataFrame, means: Dict[str, float], sds: Dict[str, float]) -> pd.DataFrame:
    Z = X.copy()
    for c in Z.columns:
        mu = float(means.get(c, Z[c].mean()))
        sd = float(sds.get(c, Z[c].std(ddof=0) or 1.0))
        if sd == 0:
            sd = 1.0
        Z[c] = (Z[c].astype(float) - mu) / sd
    return Z


def _encode_with_map(series: pd.Series, mapping: Dict[str, int]) -> Optional[np.ndarray]:
    if not mapping:
        return None
    return series.fillna("__NA__").astype(str).map(mapping).values.astype(int)


def _posterior_array(idata, name: str) -> Optional[np.ndarray]:
    try:
        arr = idata.posterior[name].to_numpy()
        return arr
    except Exception:
        return None


# ------------- RFF helpers -------------

def _rff_apply(coords: np.ndarray, params: Dict[str, Any]) -> pd.DataFrame:
    """Recreate RFF features from saved parameters (omega, b, D).
    coords: (N, d) with columns in the same order used during training.
    """
    D = int(params["D"])
    omega = np.array(params["omega"], dtype=float)  # (D, d)
    b = np.array(params["b"], dtype=float)          # (D,)
    proj = coords @ omega.T + b
    phi = np.sqrt(2.0 / D) * np.cos(proj)
    cols = [f"rff:{i}" for i in range(D)]
    return pd.DataFrame(phi, columns=cols)


def _prepare_design_with_spatial(
    table: Path,
    features_in_card: List[str],
    feature_means: Dict[str, float],
    feature_sds: Dict[str, float],
    spatial_card: Dict[str, Any],
) -> Tuple[pd.DataFrame, np.ndarray]:
    """Build X in the exact order saved in the model_card (with RFF reconstruction).
    """
    if table.suffix.lower() == ".parquet":
        df = pd.read_parquet(table)
    else:
        df = pd.read_csv(table)

    all_features = list(features_in_card)
    base_cols = [c for c in all_features if not c.startswith("rff:")]
    rff_cols = [c for c in all_features if c.startswith("rff:")]

    # Validate presence
    missing = [c for c in base_cols if c not in df.columns]
    if missing:
        raise typer.BadParameter(f"Missing base feature(s) in table: {missing}")

    X_base = df[base_cols].copy()
    if "abundance" in X_base.columns:
        if (X_base["abundance"] < -1).any():
            raise typer.BadParameter("'abundance' has values < -1; log1p undefined.")
        X_base["abundance"] = np.log1p(X_base["abundance"].astype(float))

    X_base = _standardize_with_stats(X_base, feature_means, feature_sds)

    if rff_cols:
        rff = spatial_card or {}
        if not (rff.get("rff_gp") and rff.get("rff_params")):
            raise typer.BadParameter("Model used rff:* features, but rff_params missing in model_card.")
        rff_params = rff["rff_params"]
        coord_cols = rff_params.get("coord_cols") or [c for c in ["x_um","y_um","z_um"] if c in df.columns]
        if len(coord_cols) < 2:
            raise typer.BadParameter("RFF model expects at least x_um and y_um in the table.")
        coords = df[coord_cols].to_numpy(dtype=float)
        Phi = _rff_apply(coords, rff_params)
        Phi = Phi[rff_cols]  # ensure order
        X = pd.concat([X_base, Phi], axis=1)
        X = X[all_features]
    else:
        X = X_base

    return df, X.values.astype(float)


@app.command("predict")
def predict_cmd(
    table: Path = typer.Argument(..., help="Table with features (CSV or Parquet)."),
    posterior_nc: Path = typer.Argument(..., help="posterior.nc from train."),
    model_card: Path = typer.Argument(..., help="model_card.json from train."),
    outdir: Path = typer.Option("outputs/predict", help="Output folder."),
    hdi_prob: float = typer.Option(0.95, help="Credible interval for HDIs."),
    include_draws: bool = typer.Option(False, help="Also write per-draw predictions (large)."),
    progress: bool = typer.Option(True, "--progress/--no-progress", help="Show a progress bar."),
):
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    card = json.loads(Path(model_card).read_text())
    features: List[str] = card.get("features") or []
    feature_means = card.get("feature_means") or {}
    feature_sds = card.get("feature_sds") or {}
    spatial = card.get("spatial") or {}
    group_maps = card.get("group_maps") or {}
    cage_map = group_maps.get("cage") or {}
    env_map = group_maps.get("env") or {}
    treat_map = group_maps.get("treatment") or {}

    df, X_mat = _prepare_design_with_spatial(
        table, features, feature_means, feature_sds, spatial
    )

    cage_idx = _encode_with_map(df["cage"], cage_map) if "cage" in df.columns and cage_map else None
    env_idx  = _encode_with_map(df["env_id"], env_map) if "env_id" in df.columns and env_map else None
    tr_idx   = _encode_with_map(df["treatment"], treat_map) if "treatment" in df.columns and treat_map else None

    # Load posterior
    try:
        import arviz as az
    except Exception as e:
        raise typer.BadParameter(f"arviz is required: {e}")
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
    if X_mat.shape[1] != n_feat:
        raise typer.BadParameter(f"Feature mismatch: X has {X_mat.shape[1]}, beta has {n_feat}.")

    use_progress = progress
    prog = Progress(
        SpinnerColumn(),
        TextColumn("[bold]Predicting[/bold]"),
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
        X_chunk = X_mat[start:end, :]
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

    csv_path = outdir / "predictions.csv"
    out_df.to_csv(csv_path, index=False)

    # Optional per-draw matrix (rows=observations, cols=draws)
    if include_draws:
        draws_mat = flat.T  # (N, C*D)
        npz_path = outdir / "predictions_draws.npz"
        np.savez_compressed(npz_path, draws=draws_mat)

    typer.echo(f"Predictions written to {csv_path}")
    if include_draws:
        typer.echo(f"Per-draw matrix: {npz_path}")


if __name__ == "__main__":
    app()
