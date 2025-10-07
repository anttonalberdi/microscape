from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Tuple, Optional

import yaml

# -------------------------
# Dataclasses for rule sets
# -------------------------

@dataclass
class UptakeConfig:
    mode: str = "linear_clip"       # future: "michaelis_menten"
    mM_per_uptake: float = 1.0      # 1 mM -> 1 mmol/gDW/h (magnitude)
    max_uptake: float = 10.0        # cap |lb| at this value (negative LB)
    secretion_upper: float = 1000.0 # UB for secretion

@dataclass
class TranscriptionalConfig:
    # Simple “shut lowly expressed reactions” toy config
    # Map gene IDs -> (threshold_TPM, action)
    thresholds: Dict[str, Tuple[float, str]] = field(default_factory=dict)
    # Default action if gene under threshold is encountered in a GPR
    default_action: str = "close"  # or "tighten"

@dataclass
class MetabolismRules:
    solver: str = "glpk"
    objective: Optional[str] = None
    metabolite_map: Dict[str, str] = field(default_factory=dict)
    uptake: UptakeConfig = field(default_factory=UptakeConfig)
    transcriptional: TranscriptionalConfig = field(default_factory=TranscriptionalConfig)

    # --- Environmental mapping: concentration (mM) -> exchange bounds (lb, ub) ---
    def uptake_to_bounds(self, conc_mM: float) -> Tuple[float, float]:
        """
        Convert an extracellular concentration (mM) to exchange bounds.
        Convention: uptake is negative flux (lower bound), secretion is positive (upper bound).
        """
        if self.uptake.mode == "linear_clip":
            # e.g., 5 mM -> -5 mmol/gDW/h (capped at -max_uptake)
            lb = -min(abs(conc_mM * self.uptake.mM_per_uptake), self.uptake.max_uptake)
            ub = float(self.uptake.secretion_upper)
            return lb, ub
        # Future: add Michaelis-Menten, etc.
        # Fallback
        lb = -min(abs(conc_mM * self.uptake.mM_per_uptake), self.uptake.max_uptake)
        ub = float(self.uptake.secretion_upper)
        return lb, ub

# -------------------------
# YAML I/O
# -------------------------

def load_rules(yaml_path: Path) -> MetabolismRules:
    """
    Load metabolism rules YAML (metabolism.yml) into a MetabolismRules object.
    """
    data = yaml.safe_load(Path(yaml_path).read_text()) or {}
    metab = (data.get("metabolism") or {})

    solver = metab.get("solver", "glpk")
    objective = metab.get("objective")
    metabolite_map = dict(metab.get("metabolite_map", {}))

    # uptake
    up = metab.get("uptake", {}) or {}
    uptake_cfg = UptakeConfig(
        mode=up.get("mode", "linear_clip"),
        mM_per_uptake=float(up.get("mM_per_uptake", 1.0)),
        max_uptake=float(up.get("max_uptake", 10.0)),
        secretion_upper=float(up.get("secretion_upper", 1000.0)),
    )

    # transcriptional
    tr = metab.get("transcriptional", {}) or {}
    thresholds_cfg = {}
    for gid, spec in (tr.get("thresholds", {}) or {}).items():
        # allow either scalar (threshold only) or [threshold, action]
        if isinstance(spec, (int, float)):
            thresholds_cfg[str(gid)] = (float(spec), tr.get("default_action", "close"))
        elif isinstance(spec, (list, tuple)) and len(spec) >= 1:
            thr = float(spec[0])
            act = spec[1] if len(spec) >= 2 else tr.get("default_action", "close")
            thresholds_cfg[str(gid)] = (thr, str(act))
        elif isinstance(spec, dict):
            thr = float(spec.get("threshold", 0.0))
            act = str(spec.get("action", tr.get("default_action", "close")))
            thresholds_cfg[str(gid)] = (thr, act)

    transcriptional_cfg = TranscriptionalConfig(
        thresholds=thresholds_cfg,
        default_action=str(tr.get("default_action", "close")),
    )

    return MetabolismRules(
        solver=solver,
        objective=objective,
        metabolite_map=metabolite_map,
        uptake=uptake_cfg,
        transcriptional=transcriptional_cfg,
    )

# -----------------------------------------------------------
# Optional helper: transcriptional constraints on reactions
# -----------------------------------------------------------

def expr_to_exchange_bounds(expr_tpm: float, base_lb: float, base_ub: float,
                            threshold: float, action: str = "close") -> Tuple[float, float]:
    """
    Very simple toy rule for exchange reactions driven by a single gene marker:
    If expression < threshold:
        - 'close': set lb=0, keep ub
        - 'tighten': shrink bounds (e.g., halve the range)
    Otherwise: leave as base.
    """
    if expr_tpm >= threshold:
        return base_lb, base_ub
    if action == "tighten":
        return max(0.5 * base_lb, 0.0) if base_lb < 0 else base_lb, 0.5 * base_ub
    # default 'close'
    return (0.0, base_ub)
