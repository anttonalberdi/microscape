from __future__ import annotations
import numpy as np
from pathlib import Path
import csv

from ..rd.simple import diffuse_step
from ..metabolism.mp1_stub import butyrate_rate

def run_minimal(sim_steps=200, dt=5.0, grid_shape=(1,128,128), voxel_um=10.0, progress_cb=None):
    """
    Minimal biologically-inspired demo:
    - fibre patch near the center
    - proxy production of butyrate (MP1-like), then diffusion
    - mucosa is the top band of the image (y=0..4)
    If provided, progress_cb(i, sim_steps) is called each step (or periodically by caller).
    """
    H, W = grid_shape[1], grid_shape[2]
    butyrate = np.zeros((1, H, W), dtype=float)
    fibre = np.zeros_like(butyrate)

    fibre[:, H//2-5:H//2+5, W//2-5:W//2+5] = 1.0
    cazyme = np.ones_like(butyrate) * 0.5
    anaer = np.ones_like(butyrate)

    mucosa_mask = np.zeros((1, H, W), dtype=bool)
    mucosa_mask[:, 0:5, :] = True

    for i in range(sim_steps):
        prod = butyrate_rate(cazyme, fibre, anaer, k=5e-3)
        butyrate += dt * prod
        butyrate = diffuse_step(butyrate, D=400.0, dt=dt, dx=voxel_um)
        if progress_cb is not None:
            progress_cb(i + 1, sim_steps)

    return {"butyrate": butyrate, "fibre": fibre, "mucosa_mask": mucosa_mask}

def compute_summary(butyrate: np.ndarray, mucosa_mask: np.ndarray):
    a2d = butyrate[0]
    muc = mucosa_mask[0]
    global_mean = float(np.nanmean(a2d))
    global_max  = float(np.nanmax(a2d))
    scfa_mucosa = float(np.nanmean(a2d[muc]))
    H, W = a2d.shape
    distances = np.arange(H)  # pixel rows as a proxy for distance from mucosa
    band_means = np.array([np.nanmean(a2d[y:y+1, :]) for y in range(H)], dtype=float)
    return {
        "global_mean": global_mean,
        "global_max": global_max,
        "scfa_at_mucosa": scfa_mucosa,
        "profile_dist_px": distances,
        "profile_mean": band_means,
    }

def save_summary_csv(summary: dict, out_csv: str | Path):
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        w.writerow(["global_mean", summary["global_mean"]])
        w.writerow(["global_max", summary["global_max"]])
        w.writerow(["scfa_at_mucosa", summary["scfa_at_mucosa"]])