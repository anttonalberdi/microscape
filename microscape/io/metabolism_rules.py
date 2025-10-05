# microscape/io/metabolism_rules.py
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Any
import yaml

@dataclass
class UptakeRule:
    mode: str = "linear_clip"        # currently only 'linear_clip'
    mM_per_uptake: float = 1.0       # 1 mM -> 1 mmol/gDW/h
    max_uptake: float = 10.0         # cap magnitude of uptake (|lb|)
    secretion_upper: float = 1000.0  # default UB for secretion

@dataclass
class MetabolismRules:
    solver: str = "glpk"
    objective: str | None = "BIOMASS"
    metabolite_map: Dict[str, str] = field(default_factory=dict)  # Cxxxx -> EX_... id
    uptake: UptakeRule = field(default_factory=UptakeRule)

def load_rules(yaml_path: Path) -> MetabolismRules:
    """Parse examples/demo/config/metabolism.yml (your schema)."""
    data = yaml.safe_load(Path(yaml_path).read_text()) or {}
    meta = data.get("metabolism", {}) or {}

    solver = str(meta.get("solver", "glpk")).lower()
    objective = meta.get("objective")  # may be None
    metabolite_map = dict(meta.get("metabolite_map", {}) or {})

    up = meta.get("uptake", {}) or {}
    uptake = UptakeRule(
        mode=str(up.get("mode", "linear_clip")),
        mM_per_uptake=float(up.get("mM_per_uptake", 1.0)),
        max_uptake=float(up.get("max_uptake", 10.0)),
        secretion_upper=float(up.get("secretion_upper", 1000.0)),
    )
    return MetabolismRules(
        solver=solver,
        objective=objective,
        metabolite_map=metabolite_map,
        uptake=uptake,
    )
