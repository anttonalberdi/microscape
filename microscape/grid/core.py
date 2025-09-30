from __future__ import annotations
from dataclasses import dataclass
import numpy as np

@dataclass
class Grid:
    shape: tuple[int, int, int]
    voxel_um: tuple[float, float, float]

    @property
    def coords(self):
        z = np.arange(self.shape[0]) * self.voxel_um[2]
        y = np.arange(self.shape[1]) * self.voxel_um[1]
        x = np.arange(self.shape[2]) * self.voxel_um[0]
        return np.meshgrid(x, y, z, indexing="xy")
