# microscape/runner/constraints.py
from __future__ import annotations
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
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
    lb = base_lb; ub = base_ub
    changed_env = False; changed_tx = False
    if env_lb is not None: lb = env_lb; changed_env = True
    if env_ub is not None: ub = env_ub; changed_env = True
    if tx_lb  is not None: lb = max(lb, tx_lb); changed_tx = True
    if tx_ub  is not None: ub = min(ub, tx_ub); changed_tx = True
    if lb > ub:
        mid = lb if lb < math.inf else ub
        lb = ub = mid
    mode = "both" if (changed_env and changed_tx) else ("env" if changed_env else ("tx" if changed_tx else "none"))
    return lb, ub, mode

def constrain_one(
    spot_id: str,
    microbe_id: str,
    base_bounds: Dict[str, Tuple[float,float,str]],
    env_bounds: Dict[str, Tuple[Optional[float], Optional[float]]],
    tx_bounds: Dict[str, Tuple[Optional[float], Optional[float]]],
):
    changed_ex = 0; changed_in = 0
    reaction_rows: List[Dict] = []

    for rid, (lb0, ub0, rtype) in base_bounds.items():
        env_lb, env_ub = env_bounds.get(rid, (None, None))
        tx_lb,  tx_ub  = tx_bounds.get(rid,  (None, None))
        lb_fin, ub_fin, changed = _merge_bounds(lb0, ub0, env_lb, env_ub, tx_lb, tx_ub)
        if changed != "none":
            if rtype == "exchange": changed_ex += 1
            else: changed_in += 1
        reaction_rows.append({
            "spot_id": spot_id, "microbe": microbe_id,
            "react_id": rid, "type": rtype,
            "lb_env": env_lb, "ub_env": env_ub,
            "lb_tx": tx_lb, "ub_tx": tx_ub,
            "lb_final": lb_fin, "ub_final": ub_fin,
            "changed": changed, "notes": [],
        })

    summary_row = {
        "spot_id": spot_id, "microbe": microbe_id,
        "changed_ex": int(changed_ex),
        "changed_internal": int(changed_in),
        "warnings": "",
    }
    return summary_row, reaction_rows
