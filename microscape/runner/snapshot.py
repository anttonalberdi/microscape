from __future__ import annotations
from typing import Dict, Tuple
import numpy as np
from ..io.graph_config import load_graph_yaml
from ..space.graph import build_edge_index
from ..space.diffusion import build_laplacians
from ..kinetics.sbml_engine import load_models, voxel_sources_from_sbml

def run_snapshot(config_path: str) -> Tuple[Dict[str, np.ndarray], Dict]:
    cfg = load_graph_yaml(config_path)
    nodes = cfg["space"]["nodes"]; edges = cfg["space"]["edges"]
    field_names = list(cfg["metabolite_map"].keys())

    # current concentrations per field (from 'fields' in nodes)
    C = {f: np.array([nodes[i]["fields"].get(f, 0.0) for i in range(len(nodes))], float)
         for f in field_names}

    # Laplacians (graph diffusion)
    edge_index, weights, dscale_by_field = build_edge_index(nodes, edges, field_names)
    Ls = build_laplacians(len(nodes), edge_index, weights, dscale_by_field, field_names)

    # Compute per-voxel sources/sinks from SBML (expression-constrained)
    models = load_models(cfg["models"])
    S = voxel_sources_from_sbml(cfg, models, C)  # dict f -> (N,)

    # Residual: S + D*L*C - decay*C  (decay optional)
    D_map = cfg["space"].get("diffusion_um2_s", {})
    decay_map = cfg["space"].get("decay_s", {})
    R = {}
    for f in field_names:
        D = float(D_map.get(f, 0.0))
        decay = float(decay_map.get(f, 0.0))
        R[f] = S[f].copy()
        if D:
            R[f] += D * (Ls[f] @ C[f])
        if decay:
            R[f] -= decay * C[f]

    summary = {
        "n_nodes": len(nodes),
        "fields": field_names,
        "residual_norms": {f: float(np.linalg.norm(R[f], ord=2)) for f in field_names},
        "mean_conc": {f: float(np.nanmean(C[f])) for f in field_names},
    }
    return R, summary
