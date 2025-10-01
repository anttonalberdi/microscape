from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix

@dataclass
class Edge:
    i: int
    j: int
    weight: float
    D_scale: Dict[str, float] = field(default_factory=dict)  # per-field conductance scaling

def build_laplacians(
    N: int,
    edges: List[Edge],
    field_names: List[str],
) -> Dict[str, csr_matrix]:
    """
    Build one symmetric graph Laplacian per field using per-edge conductance
    (weight * optional per-field scale).

    L_ij = L_ji = -c_ij,   L_ii += c_ij,   L_jj += c_ij
    """
    w_base = np.array([e.weight for e in edges], dtype=float)
    Ls: Dict[str, csr_matrix] = {}

    for f in field_names:
        # Per-field scaling if provided
        scale = np.array([e.D_scale.get(f, 1.0) for e in edges], dtype=float)
        conduct = w_base * scale

        ii, jj, vv = [], [], []
        for k, e in enumerate(edges):
            c = conduct[k]
            if c <= 0.0:
                continue
            i, j = e.i, e.j
            # off-diagonals
            ii.extend([i, j])
            jj.extend([j, i])
            vv.extend([-c, -c])
            # diagonals
            ii.extend([i, j])
            jj.extend([i, j])
            vv.extend([+c, +c])

        Ls[f] = coo_matrix((vv, (ii, jj)), shape=(N, N)).tocsr()

    return Ls
