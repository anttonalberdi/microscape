# microscape/io/metabolism_rules.py
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional, Any

import yaml

@dataclass
class UptakeConfig:
    mM_per_uptake: float = 1.0
    max_uptake: float = 10.0
    secretion_upper: float = 1000.0

@dataclass
class MetabolismRules:
    solver: str = "glpk"
    objective: Optional[str] = None
    metabolite_map: Dict[str, str] = field(default_factory=dict)
    uptake: UptakeConfig = field(default_factory=UptakeConfig)

def load_rules(yaml_path: Path) -> MetabolismRules:
    """
    Load config/metabolism.yml (or a path set in system.config.metabolism)
    and return a structured MetabolismRules object.
    """
    p = Path(yaml_path)
    data = yaml.safe_load(p.read_text())
    if not isinstance(data, dict) or "metabolism" not in data:
        raise ValueError(f"{p} does not define a top-level 'metabolism' section")

    m: Dict[str, Any] = data["metabolism"] or {}

    solver = str(m.get("solver", "glpk"))
    objective = m.get("objective", None)

    # map of spot metabolite IDs -> model exchange reaction IDs
    metabolite_map = dict(m.get("metabolite_map", {}) or {})

    up = m.get("uptake", {}) or {}
    uptake = UptakeConfig(
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
