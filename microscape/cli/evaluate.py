
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


def _coef_summary(arr: np.ndarray, names: List[str], hdi_prob: float = 0.95) -> pd.DataFrame:
    flat = arr.reshape(-1, arr.shape[-1])
    mean = flat.mean(axis=0)
    sd = flat.std(axis=0, ddof=0)
    lo, hi = _hdi(flat, prob=hdi_prob, axis=0)
    p_gt0 = (flat > 0).mean(axis=0)
    p_lt0 = (flat < 0).mean(axis=0)
    return pd.DataFrame({
        "name": names,
        "mean": mean,
        "sd": sd,
        f"hdi{int(hdi_prob*100)}_low": lo,
        f"hdi{int(hdi_prob*100)}_high": hi,
        "prob_gt0": p_gt0,
        "prob_lt0": p_lt0,
    })


# ------------- RFF helpers -------------

def _rff_apply(coords: np.ndarray, params: Dict[str, Any]) -> pd.DataFrame:
    """
    Recreate RFF features from saved parameters (omega, b, D).
    coords: (N, d) with columns in the same order used during training.
    """
    D = int(params["D"])
    omega = np.array(params["omega"], dtype=float)  # (D, d)
    b = np.array(params["b"], dtype=float)          # (D,)
    proj = coords @ omega.T + b
    phi = np.sqrt(2.0 / D) * np.cos(proj)
    cols = [f"rff:{i}" for i in range(D)]
    return pd.DataFrame(phi, columns=cols)


# ------------- design prep -------------

def _prepare_design_with_spatial(
    table: Path,
    features_in_card: List[str],
    target_col: str,
    feature_means: Dict[str, float],
    feature_sds: Dict[str, float],
    spatial_card: Dict[str, Any],
) -> Tuple[pd.DataFrame, np.ndarray, np.ndarray, List[str]]:
    """
    Rebuild the design matrix X in the SAME column order used by training:
      - standardize base features using saved (means, sds)
      - recreate RFF features (if used) from saved params and append in place
    Returns:
      df (original rows), X (N,K), y (N,), meta_cols (for reporting)
    """
    if table.suffix.lower() == ".parquet":
        df = pd.read_parquet(table)
    else:
        df = pd.read_csv(table)

    # Separate feature names
    all_features = list(features_in_card)
    base_cols = [c for c in all_features if not c.startswith("rff:")]
    rff_cols = [c for c in all_features if c.startswith("rff:")]

    # Validate presence of base features
    for c in [target_col, *base_cols]:
        if c not in df.columns:
            raise typer.BadParameter(f"Column missing in table: {c}")

    # Build base X with required transforms
    X_base = df[base_cols].copy()
    if "abundance" in X_base.columns:
        if (X_base["abundance"] < -1).any():
            raise typer.BadParameter("'abundance' has values < -1; log1p undefined.")
        X_base["abundance"] = np.log1p(X_base["abundance"].astype(float))

    # Standardize base cols using training stats
    X_base = _standardize_with_stats(X_base, feature_means, feature_sds)

    # Recreate RFF features if the model used them
    if rff_cols:
        rff = spatial_card or {}
        if not (rff.get("rff_gp") and rff.get("rff_params")):
            raise typer.BadParameter("Model card lists rff:* features but no rff_params found.")
        rff_params = rff["rff_params"]
        coord_cols = rff_params.get("coord_cols") or [c for c in ["x_um","y_um","z_um"] if c in df.columns]
        if len(coord_cols) < 2:
            raise typer.BadParameter("RFF model expects at least x_um and y_um in the table.")
        coords = df[coord_cols].to_numpy(dtype=float)
        Phi = _rff_apply(coords, rff_params)
        # Sanity check: number and names
        if len(rff_cols) != Phi.shape[1]:
            raise typer.BadParameter(f"RFF feature count mismatch: card={len(rff_cols)} vs computed={Phi.shape[1]}")
        Phi = Phi[rff_cols]  # order to match model
        # Concatenate base + RFF in the exact order from the card
        X = pd.concat([X_base, Phi], axis=1)
        # Ensure column order = features_in_card
        X = X[all_features]
    else:
        X = X_base

    y = df[target_col].astype(float).values
    meta_cols = [c for c in ["spot_id", "microbe", "env_id", "cage", "treatment"] if c in df.columns]
    return df, X.values.astype(float), y, meta_cols


# ---------------- CLI ----------------

@app.command("evaluate")
def evaluate_cmd(
    table: Path = typer.Argument(..., help="Training table used for the fitted model (CSV or Parquet)."),
    posterior_nc: Path = typer.Argument(..., help="Path to posterior.nc saved by 'microscape train'."),
    model_card: Path = typer.Argument(..., help="Path to model_card.json saved by 'microscape train'."),
    outdir: Path = typer.Option("outputs/evaluate", help="Where to write evaluation artifacts."),
    hdi_prob: float = typer.Option(0.95, help="Credible interval probability for HDIs."),
    progress: bool = typer.Option(True, "--progress/--no-progress", help="Show a progress bar."),
):
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # Load card
    card = json.loads(Path(model_card).read_text())
    features: List[str] = card.get("features") or []
    target_col = card.get("target") or "target"
    feature_means = card.get("feature_means") or {}
    feature_sds = card.get("feature_sds") or {}
    spatial = card.get("spatial") or {}
    group_maps = card.get("group_maps") or {}
    cage_map = group_maps.get("cage") or {}
    env_map = group_maps.get("env") or {}
    treat_map = group_maps.get("treatment") or {}

    # Build design matrix exactly like training (with RFF reconstruction if used)
    df, X_mat, y, meta_cols = _prepare_design_with_spatial(
        table, features, target_col, feature_means, feature_sds, spatial
    )

    cage_idx = _encode_with_map(df["cage"], cage_map) if "cage" in df.columns and cage_map else None
    env_idx  = _encode_with_map(df["env_id"], env_map) if "env_id" in df.columns and env_map else None
    tr_idx   = _encode_with_map(df["treatment"], treat_map) if "treatment" in df.columns and treat_map else None

    # Load posterior
    try:
        import arviz as az
    except Exception as e:
        raise typer.BadParameter(f"arviz is required to read posterior: {e}")
    idata = az.from_netcdf(posterior_nc)

    alpha = _posterior_array(idata, "alpha")
    beta  = _posterior_array(idata, "beta")
    sigma = _posterior_array(idata, "sigma")
    a_cage = _posterior_array(idata, "a_cage")
    a_env  = _posterior_array(idata, "a_env")
    gamma_t = _posterior_array(idata, "gamma_t")

    if alpha is None or beta is None:
        raise typer.BadParameter("Posterior file missing 'alpha' or 'beta'. Was the model trained correctly?")

    n_obs = df.shape[0]
    n_feat = beta.shape[-1]
    if n_feat != X_mat.shape[1]:
        raise typer.BadParameter(f"Feature mismatch: posterior has K={n_feat}, table has {X_mat.shape[1]}.")

    # Progress
    use_progress = progress
    prog = Progress(
        SpinnerColumn(),
        TextColumn("[bold]Evaluating[/bold]"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        MofNCompleteColumn(),
        TimeRemainingColumn(),
        transient=True,
    ) if use_progress else None

    if use_progress and prog is not None:
        prog.start()
        task = prog.add_task("Computing posterior means", total=n_obs)

    # Compute posterior means in chunks
    chunk = max(1, n_obs // 20)  # ~20 chunks
    mu_chunks = []
    for start in range(0, n_obs, chunk):
        end = min(n_obs, start + chunk)
        X_chunk = X_mat[start:end, :]  # (m, K)
        # alpha: (C,D), beta: (C,D,K)
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

    mu_draws = np.concatenate(mu_chunks, axis=-1)  # (C,D,N)
    mu_flat = mu_draws.reshape(-1, n_obs)
    y_hat_mean = mu_flat.mean(axis=0)
    lo, hi = _hdi(mu_flat, prob=hdi_prob, axis=0)
    residual = y - y_hat_mean

    rmse = float(np.sqrt(np.mean((y - y_hat_mean) ** 2)))
    mae = float(np.mean(np.abs(y - y_hat_mean)))
    ss_res = float(np.sum((y - y_hat_mean) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2)) or 1.0
    r2 = float(1.0 - ss_res / ss_tot)
    coverage = float(np.mean((y >= lo) & (y <= hi)))

    # Summaries
    coef_df = _coef_summary(beta, names=features, hdi_prob=hdi_prob)
    # alpha row at top
    a_flat = alpha.reshape(-1, 1)
    a_mean = a_flat.mean(axis=0)[0]
    a_sd = a_flat.std(axis=0, ddof=0)[0]
    a_lo, a_hi = _hdi(a_flat, prob=hdi_prob, axis=0)
    a_pgt0 = float((a_flat > 0).mean())
    a_plt0 = float((a_flat < 0).mean())
    coef_df.loc[-1] = ["alpha", a_mean, a_sd, a_lo[0], a_hi[0], a_pgt0, a_plt0]
    coef_df = coef_df.reset_index(drop=True)

    # Outputs
    outputs = []
    meta_cols = [c for c in ["spot_id", "microbe", "env_id", "cage", "treatment"] if c in df.columns]
    pred_df = df[meta_cols].copy() if meta_cols else pd.DataFrame(index=df.index)
    pred_df["target"] = y
    pred_df["y_hat_mean"] = y_hat_mean
    pred_df[f"y_hat_hdi{int(hdi_prob*100)}_low"] = lo
    pred_df[f"y_hat_hdi{int(hdi_prob*100)}_high"] = hi
    pred_df["residual"] = residual
    (outdir / "residuals.csv").write_text(pred_df.to_csv(index=False))
    outputs.append(("Residuals", outdir / "residuals.csv"))

    (outdir / "coef_summary.csv").write_text(coef_df.to_csv(index=False))
    outputs.append(("Coefficients", outdir / "coef_summary.csv"))

    # Metrics
    metrics = {
        "n_obs": int(n_obs),
        "hdi_prob": hdi_prob,
        "rmse": float(rmse),
        "mae": float(mae),
        "r2": float(r2),
        "coverage": float(coverage),
    }
    (outdir / "metrics.json").write_text(json.dumps(metrics, indent=2))

    typer.echo("Evaluation complete.")
    for label, path in outputs:
        typer.echo(f"  {label:18}: {path}")
    typer.echo(f"  Metrics           : {outdir / 'metrics.json'}")

if __name__ == "__main__":
    app()
