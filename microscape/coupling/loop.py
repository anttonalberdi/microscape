from __future__ import annotations
import numpy as np
from pathlib import Path
import csv

from ..rd.simple import diffuse_step
from ..metabolism.mp1_stub import butyrate_rate

def run_minimal(sim_steps=200, dt=5.0, grid_shape=(1,128,128), voxel_um=10.0):
    """
    Minimal biologically-inspired demo:
    - A fibre patch near the center
    - A simple CAZyme score (uniform)
    - Butyrate produced locally (proxy for MP1) and diffusing
    - Mucosa is the top band of the image (y=0..4)
    Returns arrays for butyrate and fibre.
    """
    H, W = grid_shape[1], grid_shape[2]
    butyrate = np.zeros((1, H, W), dtype=float)
    fibre = np.zeros_like(butyrate)

    # fibre square (10x10 px) in the middle
    fibre[:, H//2-5:H//2+5, W//2-5:W//2+5] = 1.0

    # simple scores
    cazyme = np.ones_like(butyrate) * 0.5
    anaer = np.ones_like(butyrate)       # 1.0 everywhere for demo

    # mucosa mask (top band of 5 px)
    mucosa_mask = np.zeros((1, H, W), dtype=bool)
    mucosa_mask[:, 0:5, :] = True

    for _ in range(sim_steps):
        prod = butyrate_rate(cazyme, fibre, anaer, k=5e-3)
        butyrate += dt * prod
        butyrate = diffuse_step(butyrate, D=400.0, dt=dt, dx=voxel_um)

    return {"butyrate": butyrate, "fibre": fibre, "mucosa_mask": mucosa_mask}

def compute_summary(butyrate: np.ndarray, mucosa_mask: np.ndarray):
    """
    Produce a tiny, user-friendly summary:
    - global mean/max
    - SCFA-at-mucosa score (mean butyrate within mucosa band)
    - a radial profile vs distance-from-mucosa (in pixels)
    """
    a2d = butyrate[0]
    muc = mucosa_mask[0]

    global_mean = float(np.nanmean(a2d))
    global_max = float(np.nanmax(a2d))
    scfa_mucosa = float(np.nanmean(a2d[muc]))

    # distance from top (mucosa band edge). Pixel y index serves as proxy.
    H, W = a2d.shape
    distances = np.arange(H)  # 0 at mucosa edge
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
