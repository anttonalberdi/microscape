# microscape/io/metabolism_rules.py
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional, Tuple, Any
import yaml

__all__ = [
    "UptakeConfig",
    "MetabolismRules",
    "load_rules",
    "expr_to_exchange_bounds",
    "spot_bounds_from_measurements",
]

@dataclass
class UptakeConfig:
    mode: str = "linear_clip"     # currently supported: linear_clip
    mM_per_uptake: float = 1.0    # mM -> mmol/gDW/h conversion (lower bound magnitude)
    max_uptake: float = 10.0      # cap |LB| at this value
    secretion_upper: float = 1000.0  # default UB for secretion

@dataclass
class MetabolismRules:
    solver: str = "glpk"
    objective: Optional[str] = None
    metabolite_map: Dict[str, str] = field(default_factory=dict)  # spot metabolite ID -> EX_rxn ID
    uptake: UptakeConfig = field(default_factory=UptakeConfig)

def _read_yaml(p: Path) -> dict:
    with Path(p).open("r", encoding="utf-8") as fh:
        return yaml.safe_load(fh) or {}

def load_rules(yml_path: Path) -> MetabolismRules:
    d = _read_yaml(Path(yml_path))
    m = (d.get("metabolism") or {})
    uptake = UptakeConfig(**(m.get("uptake") or {}))
    return MetabolismRules(
        solver = m.get("solver", "glpk"),
        objective = m.get("objective"),
        metabolite_map = dict(m.get("metabolite_map") or {}),
        uptake = uptake,
    )

def expr_to_exchange_bounds(
    rules: MetabolismRules,
    metabolite_id: str,
    concentration_mM: float,
) -> Optional[Tuple[str, float, float]]:
    """
    Map a *single* spot metabolite to its exchange bounds.
    Returns (exchange_rxn_id, lb, ub) or None if metabolite not mapped.
    - lb <= 0 (consumption allowed as negative flux)
    - ub >= 0 (secretion upper bound)
    """
    ex_id = rules.metabolite_map.get(metabolite_id)
    if not ex_id:
        return None
    conc = max(0.0, float(concentration_mM))
    if rules.uptake.mode == "linear_clip":
        lb_mag = conc * float(rules.uptake.mM_per_uptake)
        lb = -min(lb_mag, float(rules.uptake.max_uptake))
        ub = float(rules.uptake.secretion_upper)
    else:
        # Fallback: no consumption allowed, only secretion
        lb = 0.0
        ub = float(rules.uptake.secretion_upper)
    return ex_id, lb, ub

def spot_bounds_from_measurements(
    rules: MetabolismRules,
    metabolites_block: Dict[str, Any],
) -> Dict[str, Tuple[float, float]]:
    """
    Given a spot's `measurements.metabolites` block:
      { unit: "mM", values: { C0001: 5.2, ... } }
    return { "EX_*": (lb, ub), ... } using rules.
    """
    if not metabolites_block:
        return {}
    values = (metabolites_block.get("values") or {}) if isinstance(metabolites_block, dict) else {}
    bounds: Dict[str, Tuple[float, float]] = {}
    for met_id, conc in values.items():
        res = expr_to_exchange_bounds(rules, met_id, float(conc))
        if res:
            ex_id, lb, ub = res
            bounds[ex_id] = (lb, ub)
    return bounds
