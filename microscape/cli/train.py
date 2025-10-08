from __future__ import annotations
from pathlib import Path
from typing import List, Dict, Any, Optional
import json, typer
import numpy as np
import pandas as pd

app = typer.Typer(add_completion=False, no_args_is_help=True)

def _standardize(df: pd.DataFrame, cols: List[str]) -> pd.DataFrame:
    out = df.copy()
    for c in cols:
        mu = out[c].mean()
        sd = out[c].std(ddof=0) or 1.0
        out[c] = (out[c] - mu) / sd
    return out

@app.command("train")
def train_cmd(
    table: Path = typer.Argument(..., help="Path to model table (CSV or Parquet)."),
    outdir: Path = typer.Option("outputs/model_run", help="Output folder."),
    features: List[str] = typer.Option(None, "--feature", help="Feature columns (repeatable). Defaults: abundance + all 'met:*'."),
    group_cage: str = typer.Option("cage", help="Random-intercept group column for cages (can be missing)."),
    group_env: str = typer.Option("env_id", help="Random-intercept group column for environments (can be missing)."),
    fixed_treatment: str = typer.Option("treatment", help="Fixed-effect column for treatment (optional)."),
    draws: int = typer.Option(1000, help="Posterior draws."),
    tune: int = typer.Option(1000, help="Tuning steps."),
    chains: int = typer.Option(4, help="MCMC chains."),
    target_col: str = typer.Option("target", help="Target column."),
    seed: int = typer.Option(42, help="Random seed."),
):
    """Train a hierarchical Bayesian regression with random intercepts for cage and environment."""
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # Load table
    if table.suffix.lower()==".parquet":
        df = pd.read_parquet(table)
    else:
        df = pd.read_csv(table)

    # Default features
    if not features:
        features = []
        if "abundance" in df.columns:
            features.append("abundance")
        features += [c for c in df.columns if c.startswith("met:")]

    # Drop rows with missing target
    df = df.dropna(subset=[target_col])

    # Prepare design matrices
    X = df[features].copy()
    if "abundance" in X.columns:
        X["abundance"] = np.log1p(X["abundance"].astype(float))
    X = _standardize(X, list(X.columns))
    y = df[target_col].astype(float).values

    # Encodings
    def code_cat(series: pd.Series):
        cats = series.fillna("__NA__").astype(str).unique().tolist()
        idx = {c:i for i,c in enumerate(sorted(cats))}
        return series.fillna("__NA__").astype(str).map(idx).values.astype(int), idx

    cage_idx, cage_map = (None, {})
    if group_cage in df.columns:
        cage_idx, cage_map = code_cat(df[group_cage])
    env_idx, env_map = (None, {})
    if group_env in df.columns:
        env_idx, env_map = code_cat(df[group_env])
    treat_code, treat_map = (None, {})
    if fixed_treatment in df.columns:
        treat_code, treat_map = code_cat(df[fixed_treatment])

    try:
        import pymc as pm
        import arviz as az
    except Exception as e:
        typer.secho(f"PyMC/ArviZ are required: {e}", fg=typer.colors.RED)
        raise typer.Exit(1)

    rng = np.random.default_rng(seed)

    with pm.Model() as model:
        Xd = pm.MutableData("X", X.values)
        y_obs = pm.MutableData("y", y)

        alpha = pm.Normal("alpha", 0.0, 2.0)
        beta = pm.Normal("beta", 0.0, 1.0, shape=X.shape[1])
        sigma = pm.HalfNormal("sigma", 1.0)

        mu = alpha + pm.math.dot(Xd, beta)

        if cage_idx is not None:
            sigma_cage = pm.HalfNormal("sigma_cage", 1.0)
            a_cage = pm.Normal("a_cage", 0.0, sigma_cage, shape=len(cage_map))
            cage_d = pm.MutableData("cage_idx", cage_idx)
            mu = mu + a_cage[cage_d]

        if env_idx is not None:
            sigma_env = pm.HalfNormal("sigma_env", 1.0)
            a_env = pm.Normal("a_env", 0.0, sigma_env, shape=len(env_map))
            env_d = pm.MutableData("env_idx", env_idx)
            mu = mu + a_env[env_d]

        if treat_code is not None:
            k = len(treat_map)
            gamma_t = pm.Normal("gamma_t", 0.0, 1.0, shape=k)
            tr_d = pm.MutableData("treat_idx", treat_code)
            mu = mu + gamma_t[tr_d]

        y_like = pm.Normal("y_like", mu, sigma, observed=y_obs)

        idata = pm.sample(draws=draws, tune=tune, chains=chains, random_seed=seed, target_accept=0.9)
        ppc = pm.sample_posterior_predictive(idata, var_names=["y_like"])

    az_path = outdir / "posterior.nc"
    idata.to_netcdf(az_path)

    card = {
        "created_utc": __import__("datetime").datetime.utcnow().isoformat(),
        "table": str(table),
        "n_rows": int(df.shape[0]),
        "features": features,
        "target": target_col,
        "group_maps": {
            "cage": cage_map,
            "env": env_map,
            "treatment": treat_map,
        },
        "draws": draws,
        "tune": tune,
        "chains": chains,
        "az_path": str(az_path),
        "summary": az.summary(idata).to_dict(orient="index"),
    }
    (outdir / "model_card.json").write_text(json.dumps(card, indent=2))

    typer.echo(f"Model trained. Posterior saved to {az_path}")
    typer.echo(f"Model card: {outdir / 'model_card.json'}")
