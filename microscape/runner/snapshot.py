
from __future__ import annotations
from typing import Dict, Tuple
import numpy as np
from ..io.graph_config import load_graph_yaml
from ..space.graph import build_edge_index
from ..space.diffusion import build_laplacians
from ..kinetics.sbml_engine import load_models, voxel_sources_from_sbml

def run_snapshot(config_path: str) -> tuple[dict, dict]:
    cfg = load_graph_yaml(config_path)
    nodes = cfg["space"]["nodes"]; edges = cfg["space"]["edges"]
    field_names = list(cfg["metabolite_map"].keys())
    C = {f: np.array([nodes[i]["fields"].get(f, 0.0) for i in range(len(nodes))], float)
         for f in field_names}
    edge_index, weights, dscale_by_field = build_edge_index(nodes, edges, field_names)
    Ls = build_laplacians(len(nodes), edge_index, weights, dscale_by_field, field_names)
    models = load_models(cfg["models"])
    S = voxel_sources_from_sbml(cfg, models, C)
    D_map = cfg["space"].get("diffusion_um2_s", {})
    decay_map = cfg["space"].get("decay_s", {})
    R = {}
    for f in field_names:
        D = float(D_map.get(f, 0.0)); decay = float(decay_map.get(f, 0.0))
        Rf = S[f].copy()
        if D:
            Rf += D * (Ls[f] @ C[f])
        if decay:
            Rf -= decay * C[f]
        R[f] = Rf
    summary = {
        "n_nodes": len(nodes),
        "fields": field_names,
        "residual_norms": {f: float(np.linalg.norm(R[f], ord=2)) for f in field_names},
        "mean_conc": {f: float(np.nanmean(C[f])) for f in field_names},
    }
    return R, summary
