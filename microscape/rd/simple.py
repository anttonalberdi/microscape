from __future__ import annotations
import numpy as np

def diffuse_step(field: np.ndarray, D: float, dt: float, dx: float, decay: float=0.0):
    # simple 3D 6-neighbour explicit scheme (no-flux edges)
    lap = (
        -6*field
        + np.roll(field, 1, 0) + np.roll(field, -1, 0)
        + np.roll(field, 1, 1) + np.roll(field, -1, 1)
        + np.roll(field, 1, 2) + np.roll(field, -1, 2)
    ) / (dx*dx)
    out = field + dt * (D * lap - decay * field)
    return out
