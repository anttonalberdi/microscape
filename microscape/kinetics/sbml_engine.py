# microscape/kinetics/sbml_engine.py
from __future__ import annotations
from typing import Dict, Any
from pathlib import Path
import numpy as np
from ..core.registry import register

# Lazy globals so models are loaded once per run
_MODELS: Dict[str, Any] = {}
_MODEL_EX_RXNS: Dict[str, Dict[str, Any]] = {}  # guild -> {ex_id: rxn}
_FIELD_TO_EX: Dict[str, str] = {}               # e.g. "glc" -> "EX_glc__D_e"
_SOLVER = None

def _ensure_deps():
    try:
        import cobra  # noqa: F401
    except Exception as e:
        raise RuntimeError(
            "SBML engine requires 'cobra' (and libSBML). "
            "Install via conda-forge: conda install -c conda-forge cobra python-libsbml"
        ) from e

def _load_models(cfg: dict):
    """Load SBML models and cache exchange reactions."""
    global _MODELS, _MODEL_EX_RXNS, _FIELD_TO_EX, _SOLVER
    if _MODELS:  # already loaded
        return

    _ensure_deps()
    import cobra
    from cobra.io import read_sbml_model

    models_cfg = cfg.get("models") or {}
    if not models_cfg:
        raise RuntimeError("SBML engine: 'models:' mapping is required in YAML.")

    _SOLVER = (cfg.get("fba") or {}).get("solver", "glpk")  # glpk by default

    # field->exchange mapping used to pick exchange fluxes for tracked fields
    fmap = (cfg.get("metabolite_map") or {})
    if not fmap:
        raise RuntimeError("SBML engine: 'metabolite_map:' is required (e.g. glc: EX_glc__D_e)")
    _FIELD_TO_EX = dict(fmap)

    for guild, sbml_path in models_cfg.items():
        p = Path(sbml_path)
        if not p.exists():
            raise FileNotFoundError(f"SBML model not found for guild '{guild}': {p}")
        m = read_sbml_model(str(p))
        m.solver = _SOLVER
        # cache exchange reactions for speed
        ex_rxns = {}
        for r in m.exchanges:
            ex_rxns[r.id] = r
        _MODELS[guild] = m
        _MODEL_EX_RXNS[guild] = ex_rxns

def _node_guilds_for_idx(cfg: dict, idx: int) -> Dict[str, float]:
    """Return {guild: weight} for node idx; defaults to equal mix if absent."""
    nodes = (cfg.get("space") or {}).get("nodes") or []
    nd = nodes[idx]
    g = nd.get("guilds")
    if g:
        # normalize
        s = float(sum(g.values()))
        if s > 0:
            return {k: float(v)/s for k, v in g.items()}
    # fallback: if not provided, use all models equally
    weights = {k: 1.0 for k in (cfg.get("models") or {}).keys()}
    s = float(sum(weights.values()))
    return {k: v/s for k, v in weights.items()}

def _bounds_from_conc(field_val: float, base_uptake: float, max_uptake: float) -> float:
    """
    Map local concentration to an uptake bound (mmol/gDW/h) with a simple saturating function.
    You can later replace this with a calibrated mapping or transcript-informed scaling.
    """
    # smooth saturation: vmax * (x / (x + K)), with K derived from base_uptake
    K = max(base_uptake, 1e-6)
    return max_uptake * (field_val / (field_val + K))

def _compute_node_exchanges(node_idx: int, fields: Dict[str, np.ndarray], cfg: dict) -> Dict[str, float]:
    """
    For one voxel node, compute net exchange flux for each tracked field by combining
    guild-specific FBA solutions weighted by guild abundance.
    Returns dict field->flux (mmol/gDW/h positive = secretion).
    """
    import cobra

    guild_weights = _node_guilds_for_idx(cfg, node_idx)
    constraints = (cfg.get("constraints") or {})
    base_uptake = (constraints.get("uptake_k") or {})  # same keys as fields
    vmax = (constraints.get("max_uptake") or {})       # optional per field

    # aggregate flux over guilds
    agg_flux = {f: 0.0 for f in fields.keys()}

    for guild, w in guild_weights.items():
        model = _MODELS[guild]
        ex_rxns = _MODEL_EX_RXNS[guild]

        with model:
            # set objective — growth if present else total secretion of tracked sinks
            obj = (cfg.get("fba") or {}).get("objective", "biomass_or_sinks")
            if obj == "biomass_or_sinks":
                biomass = next((r for r in model.reactions if "biomass" in r.id.lower()), None)
                if biomass is not None:
                    model.objective = biomass
                else:
                    # sum secretion of known sinks (e.g., butyrate, propionate)
                    sinks = (cfg.get("fba") or {}).get("sink_exchanges", [])
                    model.objective = cobra.Objective(
                        sum(ex_rxns[s] for s in sinks if s in ex_rxns), direction="max"
                    )
            elif obj == "custom":
                # advanced: user can specify a linear expression, skip here
                pass

            # set uptake bounds from local concentrations
            for field, conc in fields.items():
                ex_id = _FIELD_TO_EX.get(field)
                if not ex_id or ex_id not in ex_rxns:
                    continue
                rxn = ex_rxns[ex_id]
                # COBRA convention: uptake is negative flux; set lower bound negative
                base = float(base_uptake.get(field, 0.05))
                vmax_f = float(vmax.get(field, 10.0))
                uptake_bound = _bounds_from_conc(float(conc[node_idx]), base, vmax_f)
                rxn.lower_bound = -uptake_bound  # allow uptake up to bound
                # leave rxn.upper_bound default (secretion allowed)

            sol = model.optimize()
            if sol.status != "optimal":
                # if infeasible, skip contribution
                continue

            # collect exchange fluxes for tracked fields
            for field in agg_flux.keys():
                ex_id = _FIELD_TO_EX.get(field)
                if not ex_id or ex_id not in ex_rxns:
                    continue
                fval = float(sol.fluxes.get(ex_id, 0.0))
                agg_flux[field] += w * fval

    return agg_flux  # mmol/gDW/h (positive secretion, negative uptake)

@register("sbml")
def step(fields: dict, cfg: dict, dt: float, alpha: float) -> dict:
    """
    Convert per-node exchange fluxes (mmol/gDW/h) to dC/dt per voxel field.
    'alpha' can bundle unit conversions (gDW per voxel, voxel volume, etc.).
    """
    _load_models(cfg)
    n = len(next(iter(fields.values())))  # number of nodes
    dC = {k: np.zeros(n, dtype=float) for k in fields.keys()}

    # Per-node compute
    for i in range(n):
        ex = _compute_node_exchanges(i, fields, cfg)  # mmol/gDW/h
        for field, v in ex.items():
            # map flux to concentration change: dC = alpha * flux
            dC[field][i] += alpha * v

    return dC
