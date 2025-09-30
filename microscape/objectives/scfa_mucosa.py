from __future__ import annotations
import numpy as np

def scfa_at_mucosa(butyrate: np.ndarray, mucosa_mask: np.ndarray, band_px: int = 5):
    band = butyrate[:, :band_px, :]
    return float(np.nanmean(band))
