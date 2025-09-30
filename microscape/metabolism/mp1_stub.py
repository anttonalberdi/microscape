from __future__ import annotations
import numpy as np

def butyrate_rate(cazyme: np.ndarray, fibre: np.ndarray, anaerobicity: np.ndarray, k: float=1e-3):
    c = np.clip(cazyme, 0, None)
    f = np.clip(fibre, 0, None)
    a = np.clip(anaerobicity, 0, 1)
    return k * c * f * a
