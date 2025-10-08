from __future__ import annotations
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import json, typer
import numpy as np
import pandas as pd

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _standardize_with_stats(X: pd.DataFrame, means: Dict[str, float], sds: Dict[str, float]) -> pd.DataFrame:
    Z = X.copy()
    for c in Z.columns:
        mu = float(means.get(c, Z[c].mean()))
        sd = float(sds.get(c, Z[c].std(ddof=0) or 1.0))
        if sd == 0:
            sd = 1.0
        Z[c] = (Z[c].astype(float) - mu) / sd
    return Z


def _compute_feature_stats(X: pd.DataFrame) -> Tuple[Dict[str, float], Dict[str, float]]:
    means = {c: float(X[c].mean()) for c in X.columns}
    sds = {c: float(X[c].std(ddof=0) or 1.0) for c in X.columns}
    for c in list(sds):
        if sds[c] == 0:
            sds[c] = 1.0
    return means, sds


def _validate_numeric(df: pd.DataFrame, cols: List[str], where: str = "features"):
    bad_cols = [c for c in cols if not np.issubdtype(df[c].dtype, np.number)]
    if bad_cols:
        raise typer.BadParameter(f"Non-numeric {where}: {bad_cols}. Cast them to numeric before training.")
    vals = df[cols].to_numpy()
    if not np.isfinite(vals).all():
        rows, cols_idx = np.where(~np.isfinite(vals))
        bad_col_names = sorted({cols[i] for i in cols_idx})
        raise typer.BadParameter(
            f"Found NaN/Inf in {where} after transforms. Offending columns: {bad_col_names}."
        )


def _code_cat(series: pd.Series):
    cats = series.fillna("__NA__").astype(str).unique().tolist()
    idx = {c: i for i, c in enumerate(sorted(cats))}
    return series.fillna("__NA__").astype(str).map(idx).values.astype(int), idx


@app.command("train")
def train_cmd(
    table: Path = typer.Argument(..., help="Path to model table (CSV or Parquet)."),
    outdir: Path = typer.Option("outputs/model_run", help="Output folder."),
    features: List[str] = typer.Option(
        None, "--feature",
        help="Feature columns (repeatable). Defaults: abundance + all 'met:*'."
    ),
    drop_constant: bool = typer.Option(
        True, "--drop-constant/--keep-constant",
        help="Drop features with zero variance (recommended for stability)."
    ),
    group_cage: str = typer.Option("cage", help="Random-intercept group column for cages (can be missing)."),
    group_env: str = typer.Option("env_id", help="Random-intercept group column for environments (can be missing)."),
    fixed_treatment: str = typer.Option("treatment", help="Fixed-effect column for treatment (optional)."),
    draws: int = typer.Option(1000, help="Posterior draws."),
    tune: int = typer.Option(1000, help="Tuning steps."),
    chains: int = typer.Option(4, help="MCMC chains."),
    cores: int = typer.Option(1, help="CPU cores for sampling (set 1 to reduce memory)."),
    init: str = typer.Option("adapt_diag", help="Init method: adapt_diag|jitter+adapt_diag|auto"),
    target_accept: float = typer.Option(0.9, help="NUTS target_accept (0.8–0.99)."),
    max_treedepth: int = typer.Option(12, help="Max tree depth for NUTS."),
    jax: bool = typer.Option(False, help="Use NumPyro NUTS via JAX (if installed) instead of PyMC sampler."),
    shrinkage: bool = typer.Option(False, help="Use hierarchical shrinkage prior on betas (HalfNormal global tau)."),
    student_t: bool = typer.Option(False, help="Use Student-T likelihood instead of Normal (robust to outliers)."),
    target_col: str = typer.Option("target", help="Target column."),
    seed: int = typer.Option(42, help="Random seed."),
):
    """
    Train a hierarchical Bayesian regression with random intercepts for cage and environment,
    and an optional fixed effect for treatment. Adds: input validation, cores/init/target_accept controls,
    optional JAX/NumPyro sampler, shrinkage prior, and Student-T likelihood.
    """
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # Load table
    if table.suffix.lower() == ".parquet":
        df = pd.read_parquet(table)
    else:
        df = pd.read_csv(table)

    if df.empty:
        raise typer.BadParameter("Input table is empty.")

    # Default features
    if not features:
        features = []
        if "abundance" in df.columns:
            features.append("abundance")
        features += [c for c in df.columns if c.startswith("met:")]
    # Deduplicate while preserving order
    seen = set()
    features = [c for c in features if not (c in seen or seen.add(c))]

    # Validate presence
    missing = [c for c in features if c not in df.columns]
    if missing:
        raise typer.BadParameter(f"Requested feature(s) not found in table: {missing}")

    if target_col not in df.columns:
        raise typer.BadParameter(f"Target column '{target_col}' not found in table.")

    # Drop rows with missing target
    df = df.dropna(subset=[target_col])
    if df.empty:
        raise typer.BadParameter("After dropping rows with missing target, no data remain.")

    # Build design X (apply same transform as train previously: log1p for abundance)
    X_raw = df[features].copy()
    if "abundance" in X_raw.columns:
        X_raw["abundance"] = np.log1p(X_raw["abundance"].astype(float))

    # Drop constant features if requested
    if drop_constant:
        keep = []
        dropped = []
        for c in X_raw.columns:
            v = X_raw[c].to_numpy()
            if np.nanstd(v, ddof=0) == 0:
                dropped.append(c)
            else:
                keep.append(c)
        if dropped:
            typer.secho(f"Dropping {len(dropped)} constant feature(s): {dropped}", fg=typer.colors.YELLOW)
        X_raw = X_raw[keep]
        features = keep
        if not features:
            raise typer.BadParameter("All features were constant; nothing to model.")

    # Validate numeric and finite BEFORE computing stats
    _validate_numeric(X_raw, list(X_raw.columns), where="features (after log1p)")
    y = df[target_col].astype(float).values
    if not np.isfinite(y).all():
        raise typer.BadParameter("Non-finite values in target after casting to float.")

    # Standardization stats (save for predict reproducibility)
    means, sds = _compute_feature_stats(X_raw)
    X = _standardize_with_stats(X_raw, means, sds)

    # Encodings for random/fixed factors
    cage_idx, cage_map = (None, {})
    if group_cage in df.columns:
        cage_idx, cage_map = _code_cat(df[group_cage])
    env_idx, env_map = (None, {})
    if group_env in df.columns:
        env_idx, env_map = _code_cat(df[group_env])
    treat_code, treat_map = (None, {})
    if fixed_treatment in df.columns:
        treat_code, treat_map = _code_cat(df[fixed_treatment])

    # Import PyMC/ArviZ
    try:
        import pymc as pm
        import arviz as az
    except Exception as e:
        typer.secho(f"PyMC/ArviZ are required: {e}", fg=typer.colors.RED)
        raise typer.Exit(1)

    rng = np.random.default_rng(seed)

    with pm.Model() as model:
        # Data containers (pm.Data is mutable in PyMC>=5)
        Xd = pm.Data("X", X.values.astype("float64"))   # use float64 to match default PyMC/aesara floatX
        y_obs = pm.Data("y", y.astype("float64"))

        # Priors
        alpha = pm.Normal("alpha", 0.0, 2.0)
        if shrinkage:
            tau = pm.HalfNormal("tau", 1.0)
            beta = pm.Normal("beta", 0.0, tau, shape=X.shape[1])
        else:
            beta = pm.Normal("beta", 0.0, 1.0, shape=X.shape[1])

        mu = alpha + pm.math.dot(Xd, beta)

        # Random intercepts
        if cage_idx is not None:
            sigma_cage = pm.HalfNormal("sigma_cage", 1.0)
            a_cage = pm.Normal("a_cage", 0.0, sigma_cage, shape=len(cage_map))
            cage_d = pm.Data("cage_idx", cage_idx)
            mu = mu + a_cage[cage_d]

        if env_idx is not None:
            sigma_env = pm.HalfNormal("sigma_env", 1.0)
            a_env = pm.Normal("a_env", 0.0, sigma_env, shape=len(env_map))
            env_d = pm.Data("env_idx", env_idx)
            mu = mu + a_env[env_d]

        if treat_code is not None:
            k = len(treat_map)
            gamma_t = pm.Normal("gamma_t", 0.0, 1.0, shape=k)
            tr_d = pm.Data("treat_idx", treat_code)
            mu = mu + gamma_t[tr_d]

        # Likelihood
        if student_t:
            nu = pm.Exponential("nu", 1/30)  # weakly informative heavy tails
            sigma = pm.HalfNormal("sigma", 1.0)
            y_like = pm.StudentT("y_like", nu=pm.math.softplus(nu)+1.0, mu=mu, sigma=sigma, observed=y_obs)
        else:
            sigma = pm.HalfNormal("sigma", 1.0)
            y_like = pm.Normal("y_like", mu, sigma, observed=y_obs)

        # Sampling
        if jax:
            try:
                from pymc.sampling_jax import sample_numpyro_nuts
            except Exception as e:
                raise typer.BadParameter(f"--jax requested but PyMC's sampling_jax is unavailable: {e}")
            idata = sample_numpyro_nuts(
                draws=draws, tune=tune, chains=chains,
                target_accept=target_accept, random_seed=seed,
                chain_method="vectorized"
            )
        else:
            step = pm.NUTS(max_treedepth=max_treedepth)
            idata = pm.sample(
                draws=draws, tune=tune, chains=chains,
                cores=max(1, cores), init=init, target_accept=target_accept,
                random_seed=seed, step=step
            )
            _ = pm.sample_posterior_predictive(idata, var_names=["y_like"])

    # Save posterior
    az_path = outdir / "posterior.nc"
    idata.to_netcdf(az_path)

    # Model card (add feature stats + options for reproducibility)
    card = {
        "created_utc": __import__("datetime").datetime.utcnow().isoformat(),
        "table": str(table),
        "n_rows": int(df.shape[0]),
        "features": features,
        "feature_means": {k: float(v) for k, v in _compute_feature_stats(X_raw)[0].items()},
        "feature_sds": {k: float(v) for k, v in _compute_feature_stats(X_raw)[1].items()},
        "standardization_note": "means/sds computed after log1p(abundance) and before standardization",
        "target": target_col,
        "group_maps": {
            "cage": cage_map,
            "env": env_map,
            "treatment": treat_map,
        },
        "draws": draws,
        "tune": tune,
        "chains": chains,
        "cores": max(1, cores),
        "init": init,
        "target_accept": target_accept,
        "max_treedepth": max_treedepth,
        "jax": bool(jax),
        "shrinkage": bool(shrinkage),
        "student_t": bool(student_t),
        "az_path": str(az_path),
    }

    try:
        import arviz as az
        summary = az.summary(idata)
        card["summary"] = summary.to_dict(orient="index")
    except Exception:
        pass

    (outdir / "model_card.json").write_text(json.dumps(card, indent=2))

    typer.echo(f"Model trained. Posterior saved to {az_path}")
    typer.echo(f"Model card: {outdir / 'model_card.json'}")
