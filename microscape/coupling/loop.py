from __future__ import annotations
import numpy as np
from ..rd.simple import diffuse_step
from ..metabolism.mp1_stub import butyrate_rate

def run_minimal(sim_steps=200, dt=5.0, grid_shape=(1,128,128), voxel_um=10.0):
    H, W = grid_shape[1], grid_shape[2]
    butyrate = np.zeros((1,H,W), dtype=float)

    fibre = np.zeros_like(butyrate); fibre[:, H//2-5:H//2+5, W//2-5:W//2+5] = 1.0
    cazyme = np.ones_like(butyrate) * 0.5
    anaer = np.ones_like(butyrate)

    for _ in range(sim_steps):
        prod = butyrate_rate(cazyme, fibre, anaer, k=5e-3)
        butyrate += dt * prod
        butyrate = diffuse_step(butyrate, D=400.0, dt=dt, dx=voxel_um)
    return {"butyrate": butyrate, "fibre": fibre}
