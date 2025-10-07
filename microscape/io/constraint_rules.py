from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional, Any
import yaml

@dataclass
class NormalizeSpec:
    mode: str = "per_microbe_max"     # per_microbe_max | per_model_max | constant
    constant: float = 1000.0          # used if mode == constant

@dataclass
class GPRActivitySpec:
    logic: str = "and_min_or_max"     # only option we implement now
    normalize: NormalizeSpec = field(default_factory=NormalizeSpec)
    threshold_TPM: float = 1.0
    floor_activity: float = 0.0
    cap_activity: float = 1.0

@dataclass
class BoundScalingSpec:
    min_irrev_ub: float = 0.0
    min_rev_mag: float = 0.0

@dataclass
class UptakeSpec:
    enabled: bool = True
    mM_per_uptake: float = 1.0
    max_uptake: float = 10.0
    secretion_upper: float = 1000.0
    metabolite_map: Dict[str, str] = field(default_factory=dict)

@dataclass
class WriteModelsSpec:
    enabled: bool = True
    dir: str = "constrained_models"

@dataclass
class ConstraintRules:
    solver: str = "glpk"
    gpr_activity: GPRActivitySpec = field(default_factory=GPRActivitySpec)
    bound_scaling: BoundScalingSpec = field(default_factory=BoundScalingSpec)
    uptake: UptakeSpec = field(default_factory=UptakeSpec)
    write_models: WriteModelsSpec = field(default_factory=WriteModelsSpec)

@dataclass
class UptakeRule:
    mode: str = "linear_clip"
    mM_per_uptake: float = 1.0
    max_uptake: float = 10.0
    secretion_upper: float = 1000.0

@dataclass
class TxRule:
    mode: str = "eflux"        # 'eflux' or 'threshold'
    threshold_TPM: float = 10.0
    eflux_norm: str = "max"    # 'max' or 'percentile'
    percentile: float = 95.0
    min_scale: float = 0.0     # floor scaling

def load_constraint_rules(p: Path) -> ConstraintRules:
    d = yaml.safe_load(p.read_text()) or {}
    m = d.get("metabolism", {})
    rules = ConstraintRules()
    rules.solver = m.get("solver", rules.solver)
    rules.objective = m.get("objective", rules.objective)
    rules.metabolite_map = m.get("metabolite_map", {}) or rules.metabolite_map
    u = m.get("uptake", {}) or {}
    rules.uptake = UptakeRule(
        mode=u.get("mode", rules.uptake.mode),
        mM_per_uptake=float(u.get("mM_per_uptake", rules.uptake.mM_per_uptake)),
        max_uptake=float(u.get("max_uptake", rules.uptake.max_uptake)),
        secretion_upper=float(u.get("secretion_upper", rules.uptake.secretion_upper)),
    )
    t = m.get("transcription", {}) or {}
    rules.transcription = TxRule(
        mode=t.get("mode", rules.transcription.mode),
        threshold_TPM=float(t.get("threshold_TPM", rules.transcription.threshold_TPM)),
        eflux_norm=t.get("eflux_norm", rules.transcription.eflux_norm),
        percentile=float(t.get("percentile", rules.transcription.percentile)),
        min_scale=float(t.get("min_scale", rules.transcription.min_scale)),
    )
    return rules

# ---------- Environmental bounds ----------

def env_bounds_from_conc(conc_mM: Dict[str, float], rules: ConstraintRules) -> Dict[str, Tuple[float, float]]:
    out: Dict[str, Tuple[float, float]] = {}
    for met_id, ex_id in rules.metabolite_map.items():
        c = float(conc_mM.get(met_id, 0.0))
        if c <= 0:
            out[ex_id] = (0.0, rules.uptake.secretion_upper)
        else:
            uptake = min(c * rules.uptake.mM_per_uptake, rules.uptake.max_uptake)
            out[ex_id] = (-abs(uptake), rules.uptake.secretion_upper)
    return out

# ---------- Transcriptional bounds via GPR ----------

def _parse_gpr(r: Reaction) -> Optional[str]:
    # cobra stores GPR in r.gene_reaction_rule (string), e.g. "(b0001 and b0002) or b0003"
    g = r.gene_reaction_rule
    return g if g and g.strip() else None

def _eval_gpr(expr: str, gene_to_score: Dict[str, float]) -> float:
    """
    Evaluate GPR to a single score in [0,1]:
      - 'AND' → min
      - 'OR'  → max
    """
    # Very small parser: split on 'or', then each clause split on 'and'
    # Normalize tokens and map to scores (missing genes → 0).
    def token_score(tok: str) -> float:
        t = tok.strip().strip("()").strip()
        return float(gene_to_score.get(t, 0.0))

    # OR of ANDs
    clauses = [c for c in expr.replace("AND", "and").replace("OR", "or").split("or")]
    ors = []
    for cl in clauses:
        ands = [token_score(t) for t in cl.split("and")]
        if not ands: ors.append(0.0)
        else:        ors.append(min(ands))   # AND → min
    return max(ors) if ors else 0.0          # OR → max

def _normalize_expr(gene_tpm: Dict[str, float], mode: str, percentile: float) -> float:
    if not gene_tpm:
        return 0.0
    vals = list(gene_tpm.values())
    if mode.lower() == "max":
        denom = max(vals)
    else:
        # percentile
        vs = sorted(vals)
        k = max(0, min(len(vs)-1, int(round((percentile/100.0)*(len(vs)-1)))))
        denom = vs[k]
    return float(denom) if denom > 0 else 0.0

def tx_bounds_from_transcripts(model: Model, gene_tpm: Dict[str, float], rules: ConstraintRules) -> Dict[str, Tuple[Optional[float], Optional[float]]]:
    out: Dict[str, Tuple[Optional[float], Optional[float]]] = {}
    if rules.transcription.mode.lower() == "threshold":
        thr = rules.transcription.threshold_TPM
        for r in model.reactions:
            gpr = _parse_gpr(r)
            if not gpr:
                continue
            score = _eval_gpr(gpr, {g: (1.0 if v >= thr else 0.0) for g, v in gene_tpm.items()})
            if score <= 0.0:
                # clamp both directions to zero (close reaction)
                out[r.id] = (0.0, 0.0)
        return out

    # EFlux-style
    denom = _normalize_expr(gene_tpm, rules.transcription.eflux_norm, rules.transcription.percentile)
    for r in model.reactions:
        gpr = _parse_gpr(r)
        if not gpr:
            continue
        raw = _eval_gpr(gpr, gene_tpm)
        s = rules.transcription.min_scale
        if denom > 0:
            s = max(s, min(1.0, raw / denom))
        # Scale the *magnitude* of bounds towards zero; we don't relax anything
        lb0, ub0 = float(r.lower_bound), float(r.upper_bound)
        lb_tx = lb0
        ub_tx = ub0
        if lb0 < 0: lb_tx = -abs(s*abs(lb0))
        if ub0 > 0: ub_tx =  abs(s*abs(ub0))
        out[r.id] = (lb_tx, ub_tx)
    return out

def guess_reaction_type(r: Reaction) -> str:
    # Robust exchange detection: prefix OR boundary species
    if r.id.startswith("EX_"):
        return "exchange"
    # fallback: if any external metabolite only on one side, treat as transport/exchange
    try:
        if any(m.id.endswith("_e") for m in r.metabolites):
            return "exchange"
    except Exception:
        pass
    return "internal"