
from __future__ import annotations
from typing import Dict, List
import numpy as np
import warnings, logging
import cobra

logging.getLogger("cobra").setLevel(logging.ERROR)

def load_models(model_paths: Dict[str, str]) -> Dict[str, cobra.Model]:
    models = {}
    for gid, p in model_paths.items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m = cobra.io.read_sbml_model(p)
        for r in m.reactions:
            if r.lower_bound > r.upper_bound:
                r.lower_bound, r.upper_bound = r.upper_bound, r.lower_bound
        m.solver = "glpk"
        models[gid] = m
    return models

def voxel_sources_from_sbml(cfg: dict, models: Dict[str, cobra.Model], C: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    nodes = cfg["space"]["nodes"]
    fmap = cfg["metabolite_map"]
    fba  = cfg.get("fba", {})
    cons = cfg.get("constraints", {})
    substrates = set(fba.get("substrates", []))
    sinks = set(fba.get("sinks", []))
    K = cons.get("uptake_k", {})
    Vmax = cons.get("max_uptake", {})
    alpha = float(cons.get("alpha_flux_to_dC", 1e-3))

    field_names = list(fmap.keys())
    N = len(nodes)
    S = {f: np.zeros(N, float) for f in field_names}

    for v, nd in enumerate(nodes):
        guilds = nd.get("guilds", {})
        expr_all = nd.get("expression", {})
        for gid, frac in guilds.items():
            if frac <= 0: continue
            if gid not in models: continue
            m = models[gid]
            _reset_bounds(m)
            expr = expr_all.get(gid, {})
            _apply_expression_caps(m, expr, floor=0.05)
            for f, ex_id in fmap.items():
                if ex_id not in m.reactions: continue
                rxn = m.reactions.get_by_id(ex_id)
                if ex_id in substrates:
                    C_loc = float(C[f][v])
                    cap = _monod(C_loc, float(Vmax.get(f, 0.0)), float(K.get(f, 1.0)))
                    rxn.lower_bound = -abs(cap); rxn.upper_bound = 0.0
                elif ex_id in sinks:
                    cap = float(Vmax.get(f, 1000.0))
                    rxn.lower_bound = 0.0; rxn.upper_bound = abs(cap)
            _set_objective(m, sinks if sinks else None)
            sol = m.optimize()
            if not sol or sol.status != "optimal":
                continue
            for f, ex_id in fmap.items():
                if ex_id in m.reactions:
                    flx = float(sol.fluxes.get(ex_id, 0.0))
                    S[f][v] += alpha * frac * flx
    return S

def _reset_bounds(model: cobra.Model):
    for r in model.reactions:
        if r.lower_bound > r.upper_bound:
            r.lower_bound, r.upper_bound = r.upper_bound, r.lower_bound

def _monod(C: float, Vmax: float, K: float) -> float:
    if Vmax <= 0 or C <= 0: return 0.0
    if K <= 0: return Vmax
    return Vmax * (C / (K + C))

def _apply_expression_caps(model: cobra.Model, expr: Dict[str, float], floor: float = 0.05):
    if not expr: return
    vals = list(expr.values())
    if len(vals) == 0: return
    ranks = {k: i for i, k in enumerate(sorted(expr, key=lambda k: expr[k]))}
    denom = max(len(vals)-1, 1)
    for rid, rank in ranks.items():
        a = max(rank/denom, floor)
        if rid in model.reactions:
            r = model.reactions.get_by_id(rid)
            _scale_bounds(r, a)

def _scale_bounds(rxn: cobra.Reaction, a: float):
    lb, ub = rxn.lower_bound, rxn.upper_bound
    if lb >= 0:
        rxn.upper_bound = ub * a
    elif ub <= 0:
        rxn.lower_bound = lb * a
    else:
        rxn.lower_bound = lb * a
        rxn.upper_bound = ub * a

def _set_objective(model: cobra.Model, sink_ids: List[str] | None):
    biomass = [r for r in model.reactions if "biomass" in r.id.lower()]
    coeffs = {}
    for r in biomass: coeffs[r] = 1.0
    if sink_ids:
        for sid in sink_ids:
            if sid in model.reactions:
                coeffs[model.reactions.get_by_id(sid)] = 1.0
    if coeffs:
        model.objective = coeffs
