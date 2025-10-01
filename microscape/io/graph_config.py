# add or update this loader; keep existing helpers intact
from __future__ import annotations
from pathlib import Path
import yaml

"""
Loads a YAML configuration file that describes your spatial graph (nodes/voxels and their properties)
"""

def load_graph_yaml(path: str | Path) -> dict:
    p = Path(path)
    cfg = yaml.safe_load(p.read_text())

    for nd in cfg.get("space", {}).get("nodes", []):
        if "fields" not in nd and "init" in nd:
            nd["fields"] = nd["init"]
        nd.setdefault("guilds", {})
        nd.setdefault("expression", {})

    cfg.setdefault("metabolite_map", {})
    cfg.setdefault("constraints", {})
    cfg.setdefault("fba", {})
    return cfg
