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


def _standardize_like_train(df: pd.DataFrame, feature_cols: List[str]) -> pd.DataFrame:
    X = df[feature_cols].copy()
    for c in feature_cols:
        mu = X[c].mean()
        sd = X[c].std(ddof=0) or 1.0
        X[c] = (X[c] - mu) / sd
    return X


def _prepare_design(table: Path, features: List[str], target_col: str, log1p_abundance: bool = True):
    if table.suffix.lower() == ".parquet":
        df = pd.read_parquet(table)
    else:
        df = pd.read_csv(table)

    for c in [target_col, *features]:
        if c not in df.columns:
            raise typer.BadParameter(f"Column missing in table: {c}")

    X = df[features].copy()
    if log1p_abundance and "abundance" in X.columns:
        X["abundance"] = np.log1p(X["abundance"].astype(float))
    X_std = _standardize_like_train(X, list(X.columns))

    y = df[target_col].astype(float).values
    meta_cols = [c for c in ["spot_id", "microbe", "env_id", "cage", "treatment"] if c in df.columns]

    return df, X_std, y, meta_cols


def _encode_with_map(series: pd.Series, mapping: Dict[str, int]) -> np.ndarray:
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


def _rand_eff_summary(arr: np.ndarray, level_names: List[str], hdi_prob: float = 0.95) -> pd.DataFrame:
    flat = arr.reshape(-1, arr.shape[-1])
    mean = flat.mean(axis=0)
    sd = flat.std(axis=0, ddof=0)
    lo, hi = _hdi(flat, prob=hdi_prob, axis=0)
    return pd.DataFrame({
        "level": level_names,
        "mean": mean,
        "sd": sd,
        f"hdi{int(hdi_prob*100)}_low": lo,
        f"hdi{int(hdi_prob*100)}_high": hi,
    })


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

    card = json.loads(Path(model_card).read_text())
    features: List[str] = card.get("features") or []
    target_col = card.get("target") or "target"
    group_maps = card.get("group_maps") or {}
    cage_map = group_maps.get("cage") or {}
    env_map = group_maps.get("env") or {}
    treat_map = group_maps.get("treatment") or {}

    df, X_std, y, meta_cols = _prepare_design(table, features, target_col, log1p_abundance=("abundance" in features))

    cage_idx = _encode_with_map(df["cage"], cage_map) if "cage" in df.columns and cage_map else None
    env_idx  = _encode_with_map(df["env_id"], env_map) if "env_id" in df.columns and env_map else None
    tr_idx   = _encode_with_map(df["treatment"], treat_map) if "treatment" in df.columns and treat_map else None

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
    sigma_cage = _posterior_array(idata, "sigma_cage")
    sigma_env  = _posterior_array(idata, "sigma_env")

    if alpha is None or beta is None:
        raise typer.BadParameter("Posterior file missing 'alpha' or 'beta'. Was the model trained correctly?")

    n_obs = df.shape[0]
    n_feat = beta.shape[-1]
    if n_feat != X_std.shape[1]:
        raise typer.BadParameter(f"Feature mismatch: posterior has K={n_feat}, table has {X_std.shape[1]}.")

    C, D = alpha.shape[0], alpha.shape[1]
    X_mat = X_std.values.astype(float)

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

    chunk = max(1, n_obs // 20)  # ~20 chunks
    mu_chunks = []
    for start in range(0, n_obs, chunk):
        end = min(n_obs, start + chunk)
        X_chunk = X_mat[start:end, :]  # (m, K)
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

    waic_res = None
    loo_res = None
    try:
        import arviz as az
        waic_res = az.waic(idata).to_dict()
    except Exception:
        waic_res = None
    try:
        import arviz as az
        loo_res = az.loo(idata).to_dict()
    except Exception:
        loo_res = None

    pred_df = df[meta_cols].copy() if meta_cols else pd.DataFrame(index=df.index)
    pred_df["target"] = y
    pred_df["y_hat_mean"] = y_hat_mean
    pred_df[f"y_hat_hdi{int(hdi_prob*100)}_low"] = lo
    pred_df[f"y_hat_hdi{int(hdi_prob*100)}_high"] = hi
    pred_df["residual"] = residual
    pred_df.to_csv(Path(outdir) / "residuals.csv", index=False)

    coef_df = _coef_summary(beta, names=features, hdi_prob=hdi_prob)
    coef_df.loc[-1] = ["alpha"] + [np.nan] * (coef_df.shape[1] - 1)
    coef_df.reset_index(drop=True, inplace=True)
    a_flat = alpha.reshape(-1, 1)
    a_mean = a_flat.mean(axis=0)[0]
    a_sd = a_flat.std(axis=0, ddof=0)[0]
    a_lo, a_hi = _hdi(a_flat, prob=hdi_prob, axis=0)
    a_pgt0 = float((a_flat > 0).mean())
    a_plt0 = float((a_flat < 0).mean())
    coef_df.iloc[0] = ["alpha", a_mean, a_sd, a_lo[0], a_hi[0], a_pgt0, a_plt0]
    coef_df.to_csv(Path(outdir) / "coef_summary.csv", index=False)

    if gamma_t is not None and treat_map:
        inv_map = {v: k for k, v in treat_map.items()}
        levels = [inv_map[i] for i in range(len(inv_map))]
        tr_flat = gamma_t.reshape(-1, gamma_t.shape[-1])
        tr_mean = tr_flat.mean(axis=0)
        tr_sd = tr_flat.std(axis=0, ddof=0)
        tr_lo, tr_hi = _hdi(tr_flat, prob=hdi_prob, axis=0)
        tr_df = pd.DataFrame({
            "treatment": levels,
            "mean": tr_mean,
            "sd": tr_sd,
            f"hdi{int(hdi_prob*100)}_low": tr_lo,
            f"hdi{int(hdi_prob*100)}_high": tr_hi,
        })
        tr_df.to_csv(Path(outdir) / "treatment_effects.csv", index=False)

    if a_cage is not None and cage_map:
        inv = {v: k for k, v in cage_map.items()}
        levels = [inv[i] for i in range(len(inv))]
        re_flat = a_cage.reshape(-1, a_cage.shape[-1])
        re_mean = re_flat.mean(axis=0)
        re_sd = re_flat.std(axis=0, ddof=0)
        re_lo, re_hi = _hdi(re_flat, prob=hdi_prob, axis=0)
        re_df = pd.DataFrame({
            "cage": levels,
            "mean": re_mean,
            "sd": re_sd,
            f"hdi{int(hdi_prob*100)}_low": re_lo,
            f"hdi{int(hdi_prob*100)}_high": re_hi,
        })
        re_df.to_csv(Path(outdir) / "random_intercepts_cage.csv", index=False)

    if a_env is not None and env_map:
        inv = {v: k for k, v in env_map.items()}
        levels = [inv[i] for i in range(len(inv))]
        re_flat = a_env.reshape(-1, a_env.shape[-1])
        re_mean = re_flat.mean(axis=0)
        re_sd = re_flat.std(axis=0, ddof=0)
        re_lo, re_hi = _hdi(re_flat, prob=hdi_prob, axis=0)
        re_df = pd.DataFrame({
            "env_id": levels,
            "mean": re_mean,
            "sd": re_sd,
            f"hdi{int(hdi_prob*100)}_low": re_lo,
            f"hdi{int(hdi_prob*100)}_high": re_hi,
        })
        re_df.to_csv(Path(outdir) / "random_intercepts_env.csv", index=False)

    metrics = {
        "n_obs": int(n_obs),
        "hdi_prob": hdi_prob,
        "rmse": rmse,
        "mae": mae,
        "r2": r2,
        "coverage": coverage,
        "waic": waic_res,
        "loo": loo_res,
    }
    (Path(outdir) / "metrics.json").write_text(json.dumps(metrics, indent=2))

    typer.echo("Evaluation complete.")
    typer.echo(f"  Residuals         : {Path(outdir) / 'residuals.csv'}")
    typer.echo(f"  Coefficients      : {Path(outdir) / 'coef_summary.csv'}")
    if gamma_t is not None and treat_map:
        typer.echo(f"  Treatment effects : {Path(outdir) / 'treatment_effects.csv'}")
    if a_cage is not None and cage_map:
        typer.echo(f"  RE (cage)         : {Path(outdir) / 'random_intercepts_cage.csv'}")
    if a_env is not None and env_map:
        typer.echo(f"  RE (env)          : {Path(outdir) / 'random_intercepts_env.csv'}")
    typer.echo(f"  Metrics           : {Path(outdir) / 'metrics.json'}")
