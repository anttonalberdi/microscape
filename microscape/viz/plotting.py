from __future__ import annotations
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def save_heatmap(arr: np.ndarray, path: str | Path, title: str = "", cmap: str = "viridis"):
    """
    Save a 2D heatmap for arr shaped (1,H,W) or (H,W).
    """
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    a2d = arr[0] if arr.ndim == 3 else arr
    plt.figure()
    plt.imshow(a2d, origin="lower", cmap=cmap)
    plt.title(title)
    plt.colorbar(label="a.u.")
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()

def save_profile(distances_px: np.ndarray, values: np.ndarray, path: str | Path, title: str = ""):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    plt.figure()
    plt.plot(distances_px, values, marker="o")
    plt.xlabel("Distance from mucosa (pixels)")
    plt.ylabel("Butyrate (mean)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()
