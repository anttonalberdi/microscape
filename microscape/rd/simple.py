# microscape/rd/simple.py
from __future__ import annotations
import numpy as np
import math

def _laplacian_neumann(field: np.ndarray, dx: float) -> np.ndarray:
    """
    3D 6-neighbour Laplacian with Neumann (no-flux) BC via edge padding.
    field: (Z, Y, X) or (1, Y, X) is fine.
    """
    f = np.pad(field, ((1,1),(1,1),(1,1)), mode="edge")
    c  = f[1:-1, 1:-1, 1:-1]
    xm = f[1:-1, 1:-1, 0:-2]
    xp = f[1:-1, 1:-1, 2:  ]
    ym = f[1:-1, 0:-2, 1:-1]
    yp = f[1:-1, 2:  , 1:-1]
    zm = f[0:-2, 1:-1, 1:-1]
    zp = f[2:  , 1:-1, 1:-1]
    lap = (xm + xp + ym + yp + zm + zp - 6.0 * c) / (dx * dx)
    return lap

def diffuse_step(field: np.ndarray, D: float, dt: float, dx: float, decay: float = 0.0) -> np.ndarray:
    """
    Stable explicit diffusion + linear decay with automatic sub-stepping.
    Neumann BC (no flux). Works for 2D slices shaped (1,H,W) or full 3D.
    """
    field = field.astype(np.float64, copy=True)

    # CFL stability for 3D 6-neighbour explicit scheme:
    # dt <= dx^2 / (6 * D)  (very conservative; keeps it stable)
    eps = 1e-12
    dt_max = (dx * dx) / (6.0 * max(D, eps))
    n_sub = max(1, math.ceil(dt / dt_max))
    sub_dt = dt / n_sub

    for _ in range(n_sub):
        lap = _laplacian_neumann(field, dx)
        field += sub_dt * (D * lap - decay * field)
        # optional: clamp small negatives from numerical noise
        field[field < 0] = 0.0

    return field
