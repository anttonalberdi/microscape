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
        raise typer.BadParameter(f"Non-numeric {where}: {bad_cols}. Cast to numeric before training.")
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


def _rff_features(coords: np.ndarray, D: int, lengthscale: float, seed: int) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """
    Random Fourier Features for an RBF kernel.
    coords: (N x d) with d=2 or 3 using columns [x_um, y_um, (z_um if present)]
    Returns: (DataFrame with rff:0..D-1, params dict)
    """
    rng = np.random.default_rng(seed)
    d = coords.shape[1]
    # ω ~ N(0, 1/ℓ^2 I)
    omega = rng.normal(0.0, 1.0 / max(lengthscale, 1e-9), size=(D, d))  # (D,d)
    b = rng.uniform(0.0, 2.0 * np.pi, size=D)                            # (D,)
    proj = coords @ omega.T + b                                          # (N,D)
    phi = np.sqrt(2.0 / D) * np.cos(proj)                                # (N,D)
    cols = [f"rff:{i}" for i in range(D)]
    params = {"D": int(D), "lengthscale": float(lengthscale),
              "omega": omega.tolist(), "b": b.tolist(),
              "coord_cols": None}  # filled by caller for clarity
    return pd.DataFrame(phi, columns=cols), params


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
        help="Drop features with zero variance (recommended)."
    ),
    group_cage: str = typer.Option("cage", help="Random-intercept group column for cages (can be missing)."),
    group_env: str = typer.Option("env_id", help="Random-intercept group column for environments (can be missing)."),
    fixed_treatment: str = typer.Option("treatment", help="Fixed-effect column for treatment (optional)."),
    draws: int = typer.Option(1000, help="Posterior draws."),
    tune: int = typer.Option(1000, help="Tuning (warmup) steps."),
    chains: int = typer.Option(4, help="MCMC chains."),
    cores: int = typer.Option(1, help="CPU cores (set 1 to reduce memory)."),
    init: str = typer.Option("adapt_diag", help="Init: adapt_diag|jitter+adapt_diag|auto"),
    target_accept: float = typer.Option(0.9, help="NUTS target_accept (0.8–0.99)."),
    jax: bool = typer.Option(False, help="Use NumPyro NUTS via JAX (if installed)."),
    shrinkage: bool = typer.Option(False, help="Hierarchical shrinkage on betas (global HalfNormal tau)."),
    student_t: bool = typer.Option(False, help="Student-T likelihood (robust to outliers)."),
    target_col: str = typer.Option("target", help="Target column."),
    seed: int = typer.Option(42, help="Random seed."),
    # --- spatial RFF-GP (optional) ---
    spatial_rff_gp: bool = typer.Option(False, help="Append Random Fourier Features of (x_um,y_um[,z_um]) to model a smooth spatial effect."),
    rff_D: int = typer.Option(64, help="Number of RFF features if --spatial-rff-gp."),
    rff_lengthscale: float = typer.Option(50.0, help="RBF length-scale (µm) for RFF features."),
):
    """
    Train a hierarchical Bayesian regression with random intercepts for cage and environment,
    and an optional fixed effect for treatment. Optionally add a smooth spatial effect via
    Random Fourier Features (RFF) of spot coordinates (x_um,y_um[,z_um]).
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
    # dedupe keep order
    seen = set()
    features = [c for c in features if not (c in seen or seen.add(c))]

    # Validate base feature columns
    missing = [c for c in features if c not in df.columns]
    if missing:
        raise typer.BadParameter(f"Requested feature(s) not found: {missing}")
    if target_col not in df.columns:
        raise typer.BadParameter(f"Target column '{target_col}' not found.")

    # Drop rows with missing target
    df = df.dropna(subset=[target_col])
    if df.empty:
        raise typer.BadParameter("After dropping rows with missing target, no data remain.")

    # Build design X (log1p on abundance)
    X_raw = df[features].copy()
    if "abundance" in X_raw.columns:
        if (X_raw["abundance"] < -1).any():
            badn = int((X_raw["abundance"] < -1).sum())
            raise typer.BadParameter(f"'abundance' has {badn} values < -1; log1p is undefined there.")
        X_raw["abundance"] = np.log1p(X_raw["abundance"].astype(float))

    # Drop constant features
    if drop_constant:
        keep, dropped = [], []
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

    # Validate numeric/finite before standardization
    _validate_numeric(X_raw, list(X_raw.columns), where="features (after log1p)")
    y = df[target_col].astype(float).values
    if not np.isfinite(y).all():
        raise typer.BadParameter("Non-finite values in target after casting to float.")

    # Standardize (save means/sds)
    means, sds = _compute_feature_stats(X_raw)
    X = _standardize_with_stats(X_raw, means, sds)

    # --- Spatial RFF-GP (optional) ---
    rff_params = None
    if spatial_rff_gp:
        coord_cols = [c for c in ["x_um", "y_um", "z_um"] if c in df.columns]
        if len(coord_cols) < 2:
            raise typer.BadParameter("--spatial-rff-gp requires at least x_um and y_um in the table.")
        coords = df[coord_cols].to_numpy(dtype=float)
        # Build RFF features (do NOT standardize them again; they are already scaled)
        Phi, rff_params = _rff_features(coords, D=int(rff_D), lengthscale=float(rff_lengthscale), seed=int(seed))
        rff_params["coord_cols"] = coord_cols
        # Append to X
        X = pd.concat([X, Phi.set_index(X.index)], axis=1)
        # Extend features list for bookkeeping (predict can tell these apart via name prefix)
        features = features + list(Phi.columns)

    # Encodings
    cage_idx, cage_map = (None, {})
    if "cage" in df.columns and df["cage"].notna().any():
        cage_idx, cage_map = _code_cat(df["cage"])
    env_idx, env_map = (None, {})
    if "env_id" in df.columns and df["env_id"].notna().any():
        env_idx, env_map = _code_cat(df["env_id"])
    treat_code, treat_map = (None, {})
    if "treatment" in df.columns and df["treatment"].notna().any():
        treat_code, treat_map = _code_cat(df["treatment"])

    # Import PyMC/ArviZ
    try:
        import pymc as pm
        import arviz as az
    except Exception as e:
        typer.secho(f"PyMC/ArviZ are required: {e}", fg=typer.colors.RED)
        raise typer.Exit(1)

    with pm.Model() as model:
        Xd = pm.Data("X", X.values.astype("float64"))
        y_obs = pm.Data("y", y.astype("float64"))

        # Priors
        alpha = pm.Normal("alpha", 0.0, 2.0)
        if True if not 'shrinkage' in locals() else shrinkage:  # keep original default behaviour
            pass
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
            nu = pm.Exponential("nu", 1/30)
            sigma = pm.HalfNormal("sigma", 1.0)
            pm.StudentT("y_like", nu=pm.math.softplus(nu)+1.0, mu=mu, sigma=sigma, observed=y_obs)
        else:
            sigma = pm.HalfNormal("sigma", 1.0)
            pm.Normal("y_like", mu, sigma, observed=y_obs)

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
            idata = pm.sample(
                draws=draws, tune=tune, chains=chains,
                cores=max(1, cores), init=init, target_accept=target_accept,
                random_seed=seed
            )
            _ = pm.sample_posterior_predictive(idata, var_names=["y_like"])

    # Save posterior
    az_path = outdir / "posterior.nc"
    idata.to_netcdf(az_path)

    # Model card (store feature stats + spatial params so predict can rebuild RFF)
    means_out, sds_out = _compute_feature_stats(X_raw)
    card = {
        "created_utc": __import__("datetime").datetime.utcnow().isoformat(),
        "table": str(table),
        "n_rows": int(df.shape[0]),
        "features": list(X.columns),  # includes rff:* if used
        "feature_means": {k: float(v) for k, v in means_out.items()},
        "feature_sds": {k: float(v) for k, v in sds_out.items()},
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
        "jax": bool(jax),
        "shrinkage": bool(shrinkage),
        "student_t": bool(student_t),
        "az_path": str(az_path),
        "spatial": {
            "rff_gp": bool(spatial_rff_gp),
            "rff_params": (rff_params or None),
        }
    }

    try:
        import arviz as az
        card["summary"] = az.summary(idata).to_dict(orient="index")
    except Exception:
        pass

    (outdir / "model_card.json").write_text(json.dumps(card, indent=2))

    typer.echo(f"Model trained. Posterior saved to {az_path}")
    typer.echo(f"Model card: {outdir / 'model_card.json'}")
