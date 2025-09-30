from __future__ import annotations
from pydantic import BaseModel
from pathlib import Path
import yaml

class SDPSchema(BaseModel):
    grid_voxel_um: tuple[float, float, float]
    has_masks: bool = True
    has_metabolites: bool = True

def validate_sdp(path: str | Path) -> SDPSchema:
    path = Path(path)
    meta = path / "meta" / "sample.yml"
    coords = path / "coords" / "grid.yml"
    if not meta.exists() or not coords.exists():
        raise FileNotFoundError("Missing meta/sample.yml or coords/grid.yml")
    with open(coords) as f:
        g = yaml.safe_load(f)
    schema = SDPSchema(grid_voxel_um=tuple(g.get("voxel_um", (10,10,10))),)
    return schema
