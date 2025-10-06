# microscape/io/metabolism_rules.py
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional, Tuple, Iterable
import yaml
import numpy as np

@dataclass
class UptakeConfig:
    mode: str = "linear_clip"     # linear_clip | linear | cap_only
    mM_per_uptake: float = 1.0
    max_uptake: float = 10.0
    secretion_upper: float = 1000.0

@dataclass
class TranscriptionConfig:
    mode: str = "eflux"           # eflux | threshold
    threshold_TPM: float = 10.0
    eflux_norm: str = "max"       # max | percentile
    percentile: float = 95.0      # used if eflux_norm == percentile
    min_scale: float = 0.0        # floor (0.0-1.0) when expression near zero

@dataclass
class MetabolismRules:
    solver: str = "glpk"
    objective: Optional[str] = None
    metabolite_map: Dict[str, str] = field(default_factory=dict)  # spot metabolite → EX id
    uptake: UptakeConfig = field(default_factory=UptakeConfig)
    transcription: TranscriptionConfig = field(default_factory=TranscriptionConfig)

def load_rules(yaml_path: Path) -> MetabolismRules:
    data = yaml.safe_load(Path(yaml_path).read_text())
    metab = (data.get("metabolism") or {})
    # environmental
    uptake = UptakeConfig(**(metab.get("uptake") or {}))
    # transcriptional
    transcription = TranscriptionConfig(**(metab.get("transcription") or {}))
    return MetabolismRules(
        solver = metab.get("solver", "glpk"),
        objective = metab.get("objective"),
        metabolite_map = (metab.get("metabolite_map") or {}),
        uptake = uptake,
        transcription = transcription,
    )

# ---------- Environmental mapping (mM -> EX bounds) ----------

def mM_to_uptake_bounds(mM: float, cfg: UptakeConfig) -> Tuple[float, float]:
    """
    Convert concentration (mM) to (lb, ub) for an EX reaction (uptake is negative).
    """
    if cfg.mode == "cap_only":
        lb = -min(mM, cfg.max_uptake)
    else:
        rate = mM * cfg.mM_per_uptake
        if cfg.mode == "linear_clip":
            rate = min(rate, cfg.max_uptake)
        lb = -rate
    ub = cfg.secretion_upper
    return (lb, ub)

# ---------- Transcriptional mapping (TPM -> reaction scaling) ----------

def _eval_gpr(reaction, gene_expr: Dict[str, float], threshold: float) -> float:
    """
    Compute an expression score for a reaction using GPR:
    - For AND: min(gene expr)
    - For OR:  max(gene expr)
    Fall back to max gene among reaction.genes if no rule string.
    """
    # If COBRApy parsed .genes/.gene_reaction_rule is available:
    # Use a simple heuristic: if 'and' in rule -> min; if 'or' -> max.
    rule = getattr(reaction, "gene_reaction_rule", "") or ""
    if rule:
        # tokenise: crude but OK for toy models (lowercase 'and'/'or')
        rule_l = rule.lower()
        genes = [g.id for g in getattr(reaction, "genes", [])]
        vals = [gene_expr.get(g, 0.0) for g in genes]
        if not vals:
            return 0.0
        if " and " in rule_l:
            return float(min(vals))
        if " or " in rule_l:
            return float(max(vals))
        return float(max(vals))
    else:
        genes = [g.id for g in getattr(reaction, "genes", [])]
        if not genes:
            return 0.0
        return float(max(gene_expr.get(g, 0.0) for g in genes))

def reaction_scale_from_expr(
    reaction,
    gene_expr: Dict[str, float],
    cfg: TranscriptionConfig
) -> float:
    """
    Return a [0,1] scale for reaction bounds based on expression and GPR.
    """
    score = _eval_gpr(reaction, gene_expr, cfg.threshold_TPM)

    if cfg.mode == "threshold":
        return 1.0 if score >= cfg.threshold_TPM else 0.0

    # eflux mode: scale by normalised expression
    genes = [g.id for g in getattr(reaction, "genes", [])]
    if not genes:
        return 1.0  # no GPR -> leave unconstrained

    gene_vals = np.array([gene_expr.get(g, 0.0) for g in genes], dtype=float)
    if cfg.eflux_norm == "percentile":
        denom = np.percentile(gene_vals, cfg.percentile) if len(gene_vals) else 1.0
    else:
        denom = float(gene_vals.max()) if len(gene_vals) else 1.0
    denom = max(denom, 1e-6)
    scale = float(score / denom)
    scale = max(cfg.min_scale, min(1.0, scale))
    return scale
