
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


def _posterior_array(idata, name: str) -> Optional[np.ndarray]:
    try:
        return idata.posterior[name].to_numpy()
    except Exception:
        return None


def _encode_with_map(series: pd.Series, mapping: Dict[str, int], policy: str) -> Optional[np.ndarray]:
    if not mapping:
        return None
    ser = series.fillna("__NA__").astype(str)
    idx = ser.map(mapping)
    if idx.isna().any():
        unknown = sorted(set(ser[idx.isna()]))
        if policy == "error":
            raise typer.BadParameter(f"Unknown group levels: {unknown}")
        idx = idx.fillna(-1)
    return idx.astype(int).values


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


def _parse_assignments(items: List[str]) -> Dict[str, float | str]:
    out: Dict[str, float | str] = {}
    for it in items or []:
        if "=" not in it:
            raise typer.BadParameter(f"Invalid --set '{it}', expected form col=value")
        col, val = it.split("=", 1)
        col = col.strip()
        val = val.strip()
        try:
            out[col] = float(val)
        except ValueError:
            out[col] = val
    return out


def _parse_deltas(items: List[str]) -> Dict[str, float]:
    out: Dict[str, float] = {}
    for it in items or []:
        if "=" not in it:
            raise typer.BadParameter(f"Invalid --delta '{it}', expected form col=+/-value")
        col, val = it.split("=", 1)
        col = col.strip()
        try:
            out[col] = float(val)
        except ValueError:
            raise typer.BadParameter(f"Delta value for '{col}' is not numeric: '{val}'")
    return out


def _apply_select(df: pd.DataFrame, selectors: List[str]) -> pd.DataFrame:
    if not selectors:
        return df
    mask = np.ones(len(df), dtype=bool)
    for sel in selectors:
        if "=" not in sel:
            raise typer.BadParameter(f"Invalid --select '{sel}', expected col=value")
        col, val = sel.split("=", 1)
        col = col.strip()
        val = val.strip()
        if col not in df.columns:
            raise typer.BadParameter(f"Select column '{col}' not found in table.")
        mask &= df[col].astype(str) == val
    return df.loc[mask]


@app.command("whatif")
def whatif_cmd(
    posterior_nc: Path = typer.Argument(..., help="Path to posterior.nc saved by 'microscape train'."),
    model_card: Path = typer.Argument(..., help="Path to model_card.json saved by 'microscape train'."),
    table: Path = typer.Option(None, help="Table (CSV/Parquet) to select baseline rows from. Defaults to training table in model_card."),
    outdir: Path = typer.Option("outputs/whatif", help="Where to write what-if results."),
    select: List[str] = typer.Option(None, "--select", help="Row selector(s) like 'spot_id=S0001' (ANDed). Repeatable."),
    limit: Optional[int] = typer.Option(50, help="Max number of selected rows to process (default 50)."),
    set_: List[str] = typer.Option(None, "--set", help="Set features or groups, e.g., 'met:C0010=8.0', 'abundance=100', 'treatment=diet_B'."),
    delta: List[str] = typer.Option(None, "--delta", help="Add deltas to numeric features, e.g., 'abundance=+10', 'met:C0010=-2.5'."),
    hdi_prob: float = typer.Option(0.95, help="Credible interval probability for HDIs."),
    new_level_policy: str = typer.Option("zero", help="For unseen group levels: 'zero' (no RE) or 'error'."),
    progress: bool = typer.Option(True, "--progress/--no-progress", help="Show a progress bar."),
):
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    card = json.loads(Path(model_card).read_text())
    features: List[str] = card.get("features") or []
    target_col = card.get("target") or "target"
    train_table_path = Path(card.get("table")) if card.get("table") else None
    group_maps = card.get("group_maps") or {}
    cage_map = group_maps.get("cage") or {}
    env_map = group_maps.get("env") or {}
    treat_map = group_maps.get("treatment") or {}
    feature_means = card.get("feature_means") or {}
    feature_sds = card.get("feature_sds") or {}

    base_table = Path(table) if table else train_table_path
    if base_table is None:
        raise typer.BadParameter("No table provided and training table path not found in model_card.json")
    if base_table.suffix.lower() == ".parquet":
        df_all = pd.read_parquet(base_table)
    else:
        df_all = pd.read_csv(base_table)

    df_sel = _apply_select(df_all, select)
    if limit is not None and len(df_sel) > limit:
        df_sel = df_sel.head(limit)
        typer.secho(f"Selected rows truncated to first {limit}. Use --limit to change.", fg=typer.colors.YELLOW)
    if df_sel.empty:
        raise typer.BadParameter("Selection returned 0 rows. Adjust --select filters.")

    set_vals = _parse_assignments(set_ or [])
    deltas = _parse_deltas(delta or [])

    try:
        import arviz as az
    except Exception as e:
        raise typer.BadParameter(f"arviz is required to read posterior: {e}")
    idata = az.from_netcdf(posterior_nc)

    alpha = _posterior_array(idata, "alpha")
    beta  = _posterior_array(idata, "beta")
    a_cage = _posterior_array(idata, "a_cage")
    a_env  = _posterior_array(idata, "a_env")
    gamma_t = _posterior_array(idata, "gamma_t")

    if alpha is None or beta is None:
        raise typer.BadParameter("Posterior file missing 'alpha' or 'beta'. Was the model trained correctly?")

    if feature_means and feature_sds:
        means = {k: float(v) for k, v in feature_means.items()}
        sds = {k: float(v) for k, v in feature_sds.items()}
    else:
        means, sds = _compute_train_feature_stats(train_table_path, features)

    def mu_draws_for(df_small: pd.DataFrame) -> np.ndarray:
        X_raw = df_small[features].copy()
        if "abundance" in X_raw.columns:
            if (X_raw["abundance"] < -1).any():
                raise typer.BadParameter("Found abundance < -1 after modification; log1p undefined.")
            X_raw["abundance"] = np.log1p(X_raw["abundance"].astype(float))
        X_std = _standardize_with_stats(X_raw, features, means, sds).values.astype(float)

        cage_idx = _encode_with_map(df_small["cage"], cage_map, new_level_policy) if "cage" in df_small.columns and cage_map else None
        env_idx  = _encode_with_map(df_small["env_id"], env_map, new_level_policy) if "env_id" in df_small.columns and env_map else None
        tr_idx   = _encode_with_map(df_small["treatment"], treat_map, new_level_policy) if "treatment" in df_small.columns and treat_map else None

        m = X_std.shape[0]
        base = alpha[..., None] + np.einsum("cdk,mk->cdm", beta, X_std)

        if a_cage is not None and cage_idx is not None:
            mask = cage_idx >= 0
            if mask.any():
                base[..., mask] = base[..., mask] + a_cage[:, :, cage_idx[mask]]
        if a_env is not None and env_idx is not None:
            mask = env_idx >= 0
            if mask.any():
                base[..., mask] = base[..., mask] + a_env[:, :, env_idx[mask]]
        if gamma_t is not None and tr_idx is not None:
            mask = tr_idx >= 0
            if mask.any():
                base[..., mask] = base[..., mask] + gamma_t[:, :, tr_idx[mask]]

        return base

    df_mod = df_sel.copy()
    for col, val in set_vals.items():
        if col not in df_mod.columns:
            raise typer.BadParameter(f"--set refers to column '{col}' not present in table.")
        df_mod[col] = val
    for col, d in deltas.items():
        if col not in df_mod.columns:
            raise typer.BadParameter(f"--delta refers to column '{col}' not present in table.")
        try:
            df_mod[col] = pd.to_numeric(df_mod[col], errors="raise") + d
        except Exception:
            raise typer.BadParameter(f"--delta column '{col}' is not numeric; cannot add delta.")

    use_progress = progress
    prog = Progress(
        SpinnerColumn(),
        TextColumn("[bold]What-if[/bold]"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        MofNCompleteColumn(),
        TimeRemainingColumn(),
        transient=True,
    ) if use_progress else None

    if use_progress and prog is not None:
        prog.start()
        task = prog.add_task("Computing", total=len(df_sel))

    mu_base = mu_draws_for(df_sel)
    mu_new  = mu_draws_for(df_mod)

    if use_progress and prog is not None:
        prog.advance(task, len(df_sel))
        prog.stop()

    C, D, N = mu_base.shape
    base_flat = mu_base.reshape(-1, N)
    new_flat  = mu_new.reshape(-1, N)
    delta_flat = new_flat - base_flat

    base_mean = base_flat.mean(axis=0)
    base_lo, base_hi = _hdi(base_flat, prob=hdi_prob, axis=0)

    new_mean = new_flat.mean(axis=0)
    new_lo, new_hi = _hdi(new_flat, prob=hdi_prob, axis=0)

    d_mean = delta_flat.mean(axis=0)
    d_lo, d_hi = _hdi(delta_flat, prob=hdi_prob, axis=0)

    keep_cols = [c for c in ["spot_id", "microbe", "env_id", "cage", "treatment"] if c in df_sel.columns]
    out = df_sel[keep_cols].copy() if keep_cols else pd.DataFrame(index=df_sel.index)
    out["y_base_mean"] = base_mean
    out[f"y_base_hdi{int(hdi_prob*100)}_low"] = base_lo
    out[f"y_base_hdi{int(hdi_prob*100)}_high"] = base_hi

    out["y_new_mean"] = new_mean
    out[f"y_new_hdi{int(hdi_prob*100)}_low"] = new_lo
    out[f"y_new_hdi{int(hdi_prob*100)}_high"] = new_hi

    out["delta_mean"] = d_mean
    out[f"delta_hdi{int(hdi_prob*100)}_low"] = d_lo
    out[f"delta_hdi{int(hdi_prob*100)}_high"] = d_hi

    changed_cols = sorted(set(list(_parse_assignments(set_ or {}).keys()) + list(_parse_deltas(delta or {}).keys())))
    for col in changed_cols:
        if col in df_sel.columns:
            out[f"{col}_base"] = df_sel[col].values
            out[f"{col}_new"] = df_mod[col].values

    out_csv = Path(outdir) / "whatif.csv"
    out.to_csv(out_csv, index=False)

    meta = {
        "hdi_prob": hdi_prob,
        "n_rows": int(len(out)),
        "selectors": select or [],
        "limit": limit,
        "set": set_ or [],
        "delta": delta or [],
        "new_level_policy": new_level_policy,
        "features": features,
        "uses_training_table": str(base_table == train_table_path),
        "training_table": str(train_table_path),
        "baseline_table": str(base_table),
    }
    (Path(outdir) / "whatif_meta.json").write_text(json.dumps(meta, indent=2))

    typer.echo(f"What-if results written to {out_csv}")
    typer.echo(f"Meta: {Path(outdir) / 'whatif_meta.json'}")
