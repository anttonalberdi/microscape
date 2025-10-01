
from __future__ import annotations
from pathlib import Path
from typing import Optional
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
def scatter_field(nodes_pos: np.ndarray, values: np.ndarray, out_png: Path, title: str = "", edges: Optional[np.ndarray] = None):
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(5,4), dpi=150)
    x, y = nodes_pos[:,0], nodes_pos[:,1]
    plt.scatter(x, y, c=values, s=40)
    if edges is not None:
        for i,j in edges:
            plt.plot([x[i],x[j]], [y[i],y[j]], alpha=0.2, lw=0.6, color="k")
    plt.title(title); plt.colorbar(); plt.tight_layout(); plt.savefig(out_png); plt.close()
def interpolate_to_grid(nodes_pos: np.ndarray, values: np.ndarray, out_png: Path,
                        grid_res=(128,128), method: str="linear", title: str=""):
    out_png.parent.mkdir(parents=True, exist_ok=True)
    x, y = nodes_pos[:,0], nodes_pos[:,1]
    xmin, xmax = x.min(), x.max(); ymin, ymax = y.min(), y.max()
    gx = np.linspace(xmin, xmax, grid_res[0]); gy = np.linspace(ymin, ymax, grid_res[1])
    GX, GY = np.meshgrid(gx, gy)
    Z = griddata((x, y), values, (GX, GY), method=method, fill_value=np.nan)
    plt.figure(figsize=(5,4), dpi=150)
    plt.imshow(Z, origin="lower", extent=[xmin,xmax,ymin,ymax], aspect="auto")
    plt.title(title); plt.colorbar(); plt.tight_layout(); plt.savefig(out_png); plt.close()
