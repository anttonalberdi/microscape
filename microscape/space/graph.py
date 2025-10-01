
from __future__ import annotations
from typing import Dict, List, Tuple
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix

class EdgeLike:
    def __init__(self, i, j, weight, D_scale=None):
        self.i=i; self.j=j; self.weight=weight; self.D_scale=D_scale or {}

def build_laplacians(N: int, edges, diffusion_um2_s: Dict[str, float], field_names: List[str]) -> Dict[str, csr_matrix]:
    Ls: Dict[str, csr_matrix] = {}
    w_base = np.array([e.weight for e in edges], dtype=float)
    for f in field_names:
        if any((e.D_scale or {}).get(f) is not None for e in edges):
            scale = np.array([ (e.D_scale or {}).get(f, 1.0) for e in edges ], dtype=float)
        else:
            scale = 1.0
        conduct = w_base * scale
        ii, jj, vv = [], [], []
        for k, e in enumerate(edges):
            i, j = e.i, e.j
            c = conduct[k]
            if c == 0: continue
            ii += [i, j, i, j]
            jj += [j, i, i, j]
            vv += [-c, -c, c, c]
        Ls[f] = coo_matrix((vv, (ii, jj)), shape=(N, N)).tocsr()
    return Ls

def explicit_diffusion_step(fields: Dict[str, np.ndarray], Ls: Dict[str, csr_matrix],
                            D_map: Dict[str, float], dt: float):
    for f, arr in fields.items():
        D = float(D_map.get(f, 0.0))
        if D:
            fields[f] = arr + dt * (D * (Ls[f] @ arr))

def stable_dt_upper_bound(L: csr_matrix, D: float, safety: float = 0.45) -> float:
    if D <= 0: return 1e9
    abs_rowsum = np.abs(L).sum(axis=1).A.ravel().max()
    if abs_rowsum == 0: return 1e9
    return safety / (D * abs_rowsum)

def coarsen(nodes, edges, groups: List[List[int]]):
    mapping = {}
    for new_i, grp in enumerate(groups):
        for old_i in grp: mapping[old_i] = new_i
    new_nodes = []
    for new_i, grp in enumerate(groups):
        vol = sum(nodes[i].volume_nl for i in grp)
        rep = nodes[grp[0]]
        new_nodes.append(type(rep)(id=f"g{new_i}", pos=rep.pos, volume_nl=vol,
                                   init={}, guilds=rep.guilds, transcripts=rep.transcripts, region=rep.region))
    E = {}
    for e in edges:
        i2, j2 = mapping[e.i], mapping[e.j]
        if i2 == j2: continue
        key = tuple(sorted((i2, j2)))
        if key not in E:
            E[key] = {"i": key[0], "j": key[1], "weight": 0.0, "D_scale": {}}
        E[key]["weight"] += e.weight
        for k,v in (e.D_scale or {}).items():
            E[key]["D_scale"][k] = E[key]["D_scale"].get(k, 0.0) + float(v)
    new_edges = [EdgeLike(**rec) for rec in E.values()]
    return new_nodes, new_edges
