from __future__ import annotations
from typing import Dict
import numpy as np
from scipy.sparse import csr_matrix

def explicit_diffusion_step(
    fields: Dict[str, np.ndarray],
    Ls: Dict[str, csr_matrix],
    D_map: Dict[str, float],
    dt: float,
) -> Dict[str, np.ndarray]:
    """
    u(t+dt) = u + dt * ( D * L * u )
    Returns a NEW dict to avoid in-place surprises.
    """
    out = {}
    for f, u in fields.items():
        D = float(D_map.get(f, 0.0))
        if D == 0.0:
            out[f] = u
        else:
            out[f] = u + dt * (D * (Ls[f] @ u))
    return out

def stable_dt_upper_bound(L: csr_matrix, D: float, safety: float = 0.45) -> float:
    """
    Crude but conservative bound based on max |row sum| for explicit Euler stability.
    dt <= safety / (D * max_row_sum)
    """
    if D <= 0:
        return 1e9
    abs_rowsum = np.abs(L).sum(axis=1).A.ravel().max()
    if abs_rowsum == 0:
        return 1e9
    return safety / (D * abs_rowsum)
