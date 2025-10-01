from __future__ import annotations
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix

def build_laplacians(N: int, edge_index, weights, dscale_by_field, field_names):
    Ls = {}
    i = edge_index[0]; j = edge_index[1]
    for f in field_names:
        w = weights * dscale_by_field[f]
        # off-diagonal entries
        data = -w
        # degree (row-sum of conductances)
        deg = np.bincount(i, w, minlength=N)
        # assemble
        L = coo_matrix((data, (i, j)), shape=(N, N)).tocsr()
        L.setdiag(L.diagonal() + deg)
        Ls[f] = L
    return Ls
