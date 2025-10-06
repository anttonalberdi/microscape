# microscape/runner/constraints.py
from __future__ import annotations
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from pathlib import Path
import math

@dataclass
class ReactionBound:
    react_id: str
    rtype: str             # 'exchange' | 'internal'
    lb_env: Optional[float]
    ub_env: Optional[float]
    lb_tx: Optional[float]
    ub_tx: Optional[float]
    lb_final: float
    ub_final: float
    changed: str           # 'env'|'tx'|'both'|'none'
    notes: List[str]

def _merge_bounds(
    base_lb: float, base_ub: float,
    env_lb: Optional[float], env_ub: Optional[float],
    tx_lb: Optional[float], tx_ub: Optional[float],
) -> Tuple[float, float, str]:
    lb = base_lb
    ub = base_ub
    changed_env = False
    changed_tx = False
    if env_lb is not None: lb, changed_env = env_lb, True
    if env_ub is not None: ub, changed_env = env_ub, True if env_lb is None else True
    if tx_lb  is not None: lb, changed_tx  = max(lb, tx_lb), True   # tighten (don’t relax)
    if tx_ub  is not None: ub, changed_tx  = min(ub, tx_ub), True
    # Ensure feasibility ordering
    if lb > ub:
        # fallback: clamp to equality at the tighter one
        mid = lb if lb < math.inf else ub
        lb = ub = mid
    if changed_env and changed_tx: mode = "both"
    elif changed_env: mode = "env"
    elif changed_tx: mode = "tx"
    else: mode = "none"
    return lb, ub, mode

def constrain_one(
    spot_id: str,
    microbe_id: str,
    # base model info for reaction types and default bounds:
    base_bounds: Dict[str, Tuple[float,float,str]],  # react_id -> (lb_default, ub_default, 'exchange'|'internal')
    # proposed constraints from each source:
    env_bounds: Dict[str, Tuple[Optional[float], Optional[float]]],  # react_id -> (lb_env, ub_env)
    tx_bounds: Dict[str, Tuple[Optional[float], Optional[float]]],   # react_id -> (lb_tx,  ub_tx)
) -> Tuple[Dict, List[Dict]]:
    """
    Returns:
      summary_row: dict(spot_id, microbe, mode, changed_ex, changed_internal, warnings)
      reaction_rows: list of dict per reaction (for JSON assembly)
    """
    changed_ex = 0
    changed_in = 0
    warnings: List[str] = []
    reaction_rows: List[Dict] = []

    for rid, (lb0, ub0, rtype) in base_bounds.items():
        env_lb, env_ub = env_bounds.get(rid, (None, None))
        tx_lb,  tx_ub  = tx_bounds.get(rid,  (None, None))

        lb_fin, ub_fin, changed = _merge_bounds(lb0, ub0, env_lb, env_ub, tx_lb, tx_ub)
        if changed != "none":
            if rtype == "exchange": changed_ex += 1
            else: changed_in += 1

        reaction_rows.append({
            "spot_id": spot_id,
            "microbe": microbe_id,
            "react_id": rid,
            "type": rtype,
            "lb_env": env_lb,
            "ub_env": env_ub,
            "lb_tx": tx_lb,
            "ub_tx": tx_ub,
            "lb_final": lb_fin,
            "ub_final": ub_fin,
            "changed": changed,
            "notes": [],
        })

    summary_row = {
        "spot_id": spot_id,
        "microbe": microbe_id,
        "changed_ex": int(changed_ex),
        "changed_internal": int(changed_in),
        "warnings": ";".join(warnings) if warnings else "",
    }
    return summary_row, reaction_rows
