
from __future__ import annotations
import numpy as np

def build_edge_index(nodes, edges, field_names):
    id_to_idx = {nd["id"]: i for i, nd in enumerate(nodes)}
    ii, jj, w = [], [], []
    dscale = {f: [] for f in field_names}
    for e in edges:
        i = id_to_idx[e["i"]]; j = id_to_idx[e["j"]]
        wt = float(e.get("weight", 1.0))
        ii += [i, j]; jj += [j, i]; w += [wt, wt]
        for f in field_names:
            s = (e.get("D_scale", {}) or {}).get(f, 1.0)
            dscale[f] += [s, s]
    edge_index = np.vstack([ii, jj]).astype(int)
    weights = np.array(w, float)
    dscale_by_field = {f: np.array(dscale[f], float) for f in field_names}
    return edge_index, weights, dscale_by_field
