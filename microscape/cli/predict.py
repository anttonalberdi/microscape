
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


def _hdi(x: np.ndarray, prob: float = 0.95, axis: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
    lo = (1.0 - prob) / 2.0
    hi = 1.0 - lo
    q = np.quantile(x, [lo, hi], axis=axis)
    return q[0], q[1]


def _standardize_with_stats(df: pd.DataFrame, feature_cols: List[str], means: Dict[str, float], sds: Dict[str, float]) -> pd.DataFrame:
    X = df[feature_cols].copy()
    for c in feature_cols:
        mu = float(means.get(c, X[c].mean()))
        sd = float(sds.get(c, X[c].std(ddof=0) or 1.0))
        if sd == 0:
            sd = 1.0
        X[c] = (X[c].astype(float) - mu) / sd
    return X


def _compute_train_feature_stats(train_table: Path, features: List[str]) -> Tuple[Dict[str, float], Dict[str, float]]:
    if train_table.suffix.lower() == ".parquet":
        tr = pd.read_parquet(train_table)
    else:
        tr = pd.read_csv(train_table)
    X = tr[features].copy()
    if "abundance" in X.columns:
        X["abundance"] = np.log1p(X["abundance"].astype(float))
    means = {c: float(X[c].mean()) for c in X.columns}
    sds = {c: float(X[c].std(ddof=0) or 1.0) for c in X.columns}
    for c in list(sds):
        if sds[c] == 0:
            sds[c] = 1.0
    return means, sds


def _encode_with_map(series: pd.Series, mapping: Dict[str, int], policy: str) -> Optional[np.ndarray]:
    """
    Map categorical labels to indices using a stored mapping from training.
    policy: 'zero' -> unknown levels get None (contributes 0 RE)
            'error'-> raise if unknown levels appear
    """
    if not mapping:
        return None
    ser = series.fillna("__NA__").astype(str)
    idx = ser.map(mapping)
    if idx.isna().any():
        unknown = sorted(set(ser[idx.isna()]))
        if policy == "error":
            raise typer.BadParameter(f"Unknown levels for grouping variable: {unknown}. Use --new-level-policy zero to treat as 0 random effect.")
        # 'zero': return array with -1 for unknowns and handle later by skipping RE
        out = idx.fillna(-1).astype(int).values
        return out
    return idx.astype(int).values


def _posterior_array(idata, name: str) -> Optional[np.ndarray]:
    try:
        arr = idata.posterior[name].to_numpy()
        return arr
    except Exception:
        return None


@app.command("predict")
def predict_cmd(
    posterior_nc: Path = typer.Argument(..., help="Path to posterior.nc saved by 'microscape train'."),
    model_card: Path = typer.Argument(..., help="Path to model_card.json saved by 'microscape train'."),
    table: Path = typer.Option(None, help="Table (CSV/Parquet) with features to predict on. Defaults to training table in model_card."),
    outdir: Path = typer.Option("outputs/predict", help="Where to write predictions."),
    hdi_prob: float = typer.Option(0.95, help="Credible interval probability for prediction HDIs (from mu)."),
    ppc: bool = typer.Option(False, "--ppc/--no-ppc", help="If set, also simulate predictive distribution using sigma."),
    new_level_policy: str = typer.Option("zero", help="For unseen group levels: 'zero' (no RE) or 'error'."),
    progress: bool = typer.Option(True, "--progress/--no-progress", help="Show a progress bar."),
):
    """
    Generate posterior predictions for rows in a modeling table.

    - Features are standardized using training means/SDs (inferred from training table path in model_card).
    - Group encodings (cage/env/treatment) follow training maps in model_card.
    - Unseen group levels are treated as 0 random effect by default (configurable).
    """
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # Load model card
    card = json.loads(Path(model_card).read_text())
    features: List[str] = card.get("features") or []
    target_col = card.get("target") or "target"
    train_table_path = Path(card.get("table")) if card.get("table") else None
    group_maps = card.get("group_maps") or {}
    cage_map = group_maps.get("cage") or {}
    env_map = group_maps.get("env") or {}
    treat_map = group_maps.get("treatment") or {}

    # Load training stats for standardization
    if train_table_path is None:
        raise typer.BadParameter("model_card.json is missing 'table' path to the training table; required to recover standardization stats.")
    means, sds = _compute_train_feature_stats(train_table_path, features)

    # Load table for prediction
    pred_table = Path(table) if table else train_table_path
    if pred_table.suffix.lower() == ".parquet":
        df = pd.read_parquet(pred_table)
    else:
        df = pd.read_csv(pred_table)

    # Validate features
    for c in features:
        if c not in df.columns:
            raise typer.BadParameter(f"Feature '{c}' missing from prediction table {pred_table}.")

    # Build design matrix
    X = df[features].copy()
    if "abundance" in X.columns:
        X["abundance"] = np.log1p(X["abundance"].astype(float))
    X_std = _standardize_with_stats(X, features, means, sds)
    X_mat = X_std.values.astype(float)

    # Encode groups
    cage_idx = _encode_with_map(df["cage"], cage_map, new_level_policy) if "cage" in df.columns and cage_map else None
    env_idx  = _encode_with_map(df["env_id"], env_map, new_level_policy) if "env_id" in df.columns and env_map else None
    tr_idx   = _encode_with_map(df["treatment"], treat_map, new_level_policy) if "treatment" in df.columns and treat_map else None

    # Load posterior
    try:
        import arviz as az
    except Exception as e:
        raise typer.BadParameter(f"arviz is required to read posterior: {e}")
    idata = az.from_netcdf(posterior_nc)

    alpha = _posterior_array(idata, "alpha")             # (C,D)
    beta  = _posterior_array(idata, "beta")              # (C,D,K)
    a_cage = _posterior_array(idata, "a_cage")           # (C,D,Jc) or None
    a_env  = _posterior_array(idata, "a_env")            # (C,D,Je) or None
    gamma_t = _posterior_array(idata, "gamma_t")         # (C,D,Jt) or None
    sigma = _posterior_array(idata, "sigma")             # (C,D)  for PPC

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
        TextColumn("[bold]Predicting[/bold]"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        MofNCompleteColumn(),
        TimeRemainingColumn(),
        transient=True,
    ) if use_progress else None

    if use_progress and prog is not None:
        prog.start()
        task = prog.add_task("Computing posterior means", total=n_obs)

    # Compute mu draws in chunks
    C, D = alpha.shape[0], alpha.shape[1]
    chunk = max(1, n_obs // 20)
    mu_chunks = []
    for start in range(0, n_obs, chunk):
        end = min(n_obs, start + chunk)
        X_chunk = X_mat[start:end, :]   # (m, K)
        base = alpha[..., None] + np.einsum("cdk,mk->cdm", beta, X_chunk)

        # Random effects; handle unknown levels (encoded as -1) by skipping addition
        if a_cage is not None and cage_idx is not None:
            idx = cage_idx[start:end]  # (m,)
            if (idx >= 0).any():
                # For valid indices add RE; unknown (-1) contribute 0
                mask = idx >= 0
                if mask.any():
                    # add only to rows with known levels
                    base[..., mask] = base[..., mask] + a_cage[:, :, idx[mask]]
        if a_env is not None and env_idx is not None:
            idx = env_idx[start:end]
            if (idx >= 0).any():
                mask = idx >= 0
                if mask.any():
                    base[..., mask] = base[..., mask] + a_env[:, :, idx[mask]]
        if gamma_t is not None and tr_idx is not None:
            idx = tr_idx[start:end]
            if (idx >= 0).any():
                mask = idx >= 0
                if mask.any():
                    base[..., mask] = base[..., mask] + gamma_t[:, :, idx[mask]]

        mu_chunks.append(base)
        if use_progress and prog is not None:
            prog.advance(task, end - start)

    if use_progress and prog is not None:
        prog.stop()

    mu_draws = np.concatenate(mu_chunks, axis=-1)  # (C,D,N)
    mu_flat = mu_draws.reshape(-1, n_obs)

    # Summaries (conditional means)
    y_hat_mean = mu_flat.mean(axis=0)
    lo, hi = _hdi(mu_flat, prob=hdi_prob, axis=0)

    pred_df = pd.DataFrame({
        "y_hat_mean": y_hat_mean,
        f"y_hat_hdi{int(hdi_prob*100)}_low": lo,
        f"y_hat_hdi{int(hdi_prob*100)}_high": hi,
    })

    # Include meta columns if present
    for c in ["spot_id", "microbe", "env_id", "cage", "treatment"]:
        if c in df.columns:
            pred_df[c] = df[c].values

    # Predictive distribution (optional): simulate y_like using sigma
    if ppc and sigma is not None:
        # Sample one noise draw per posterior draw for each observation would be heavy.
        # Instead, compute predictive HDIs by adding N(0, sigma) to mu per draw via Monte Carlo (moderate size).
        # We'll use 200 draws uniformly across chains.
        draws = mu_flat.shape[0]
        use = min(200, draws)
        sel = np.linspace(0, draws - 1, use, dtype=int)
        mu_sel = mu_flat[sel, :]  # (use, N)
        # sigma per draw: flatten (C*D,)-> (draws,), select same indices
        sig_flat = sigma.reshape(-1)
        sig_sel = sig_flat[sel]    # (use,)
        # noise ~ N(0, sigma)
        rng = np.random.default_rng(2025)
        noise = rng.normal(0.0, sig_sel[:, None], size=mu_sel.shape)
        y_like = mu_sel + noise
        plo, phi = _hdi(y_like, prob=hdi_prob, axis=0)
        pred_df[f"y_pred_hdi{int(hdi_prob*100)}_low"] = plo
        pred_df[f"y_pred_hdi{int(hdi_prob*100)}_high"] = phi

    # Write
    out_csv = Path(outdir) / "predictions.csv"
    pred_df.to_csv(out_csv, index=False)

    # Small metadata
    meta = {
        "hdi_prob": hdi_prob,
        "ppc": bool(ppc),
        "n_rows": int(len(pred_df)),
        "uses_training_table": str(pred_table == train_table_path),
        "training_table": str(train_table_path),
        "prediction_table": str(pred_table),
        "features": features,
        "new_level_policy": new_level_policy,
    }
    (Path(outdir) / "predict_meta.json").write_text(json.dumps(meta, indent=2))

    typer.echo(f"Predictions written to {out_csv}")
    typer.echo(f"Meta: {Path(outdir) / 'predict_meta.json'}")
