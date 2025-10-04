from __future__ import annotations
from pathlib import Path
from typing import Dict, Any
import yaml

def _max_expr_for_any(genes: Dict[str, float], candidates: list[str]) -> float:
    if not genes or not candidates:
        return 0.0
    return max((float(genes.get(g, 0.0)) for g in candidates), default=0.0)

def _score(expr: float, thresh: float, mode: str = "clip_linear") -> float:
    if thresh <= 0:
        return 0.0
    if mode == "clip_linear":
        s = expr / thresh
        return 1.0 if s > 1.0 else (0.0 if s < 0.0 else s)
    # fallback = binary
    return 1.0 if expr >= thresh else 0.0

def _validate_rules(rules: dict) -> None:
    for key in ("thresholds", "markers"):
        if key not in rules:
            raise ValueError(f"Missing '{key}' in ecology rules.")
    th = rules["thresholds"]
    mrk = rules["markers"]
    for req in ("adhesion_TPM", "motility_TPM", "stress_TPM"):
        if req not in th:
            raise ValueError(f"Missing thresholds.{req}")
    for trait in ("adhesion", "motility", "stress"):
        if trait not in mrk or "genes" not in mrk[trait]:
            raise ValueError(f"Missing markers.{trait}.genes")
        if not isinstance(mrk[trait]["genes"], list):
            raise ValueError(f"markers.{trait}.genes must be a list of gene codes")

# microscape/profile/ecology.py (core changes shown)
def infer_ecology_for_spot(spot_yaml: Path, rules_yaml: Path) -> Dict[str, Any]:
    spot = yaml.safe_load(Path(spot_yaml).read_text())
    cfg  = yaml.safe_load(Path(rules_yaml).read_text())["ecology"]

    scoring_mode = (cfg.get("scoring") or {}).get("mode", "clip_linear")
    traits = cfg["traits"]

    tblock = ((spot.get("measurements") or {}).get("transcripts") or {})
    if (tblock.get("schema") or "genes") != "genes":
        raise ValueError("transcripts.schema must be 'genes'")
    per_microbe = tblock.get("values", {})

    out = {}
    for m_id, gene_expr in per_microbe.items():
        m_states = {}
        for tr in traits:
            tr_id   = tr["id"]
            states  = tr["states"]; assert len(states) == 2
            pos     = tr.get("positive_state", states[0])
            # overrides
            ov      = (tr.get("overrides") or {}).get(m_id, {})
            markers = ov.get("markers", tr["markers"])
            thr     = float(ov.get("threshold", tr["threshold"]))

            expr = max(float(gene_expr.get(g, 0.0)) for g in markers) if markers else 0.0

            if scoring_mode == "clip_linear":
                score = max(0.0, min(1.0, expr / thr if thr > 0 else 0.0))
            elif scoring_mode == "binary":
                score = 1.0 if expr >= thr else 0.0
            elif scoring_mode == "logistic":
                pars = (cfg.get("scoring") or {}).get("logistic", {})
                k = float(pars.get("k", 1.0))
                # center at threshold: logistic((expr-thr))
                import math
                score = 1.0 / (1.0 + math.exp(-k * (expr - thr)))
            else:
                raise ValueError(f"Unknown scoring mode: {scoring_mode}")

            state = pos if score >= 0.5 else (states[1] if pos == states[0] else states[0])
            m_states[tr_id] = {
                "state": state,
                "score": round(score, 3),
                "marker_max_TPM": round(expr, 3),
                "markers": markers,
                "threshold": thr,
            }
        out[m_id] = m_states
    return out

def write_enriched_spot(spot_yaml: Path, rules_yaml: Path, out_yaml: Path) -> None:
    spot = yaml.safe_load(Path(spot_yaml).read_text())
    ecology = infer_ecology_for_spot(spot_yaml, rules_yaml)
    spot.setdefault("inference", {})
    spot["inference"]["ecology"] = {
        "method": "rules-v1-codes",
        "rules_file": str(rules_yaml),
        "states": ecology,
    }
    Path(out_yaml).parent.mkdir(parents=True, exist_ok=True)
    Path(out_yaml).write_text(yaml.safe_dump(spot, sort_keys=False))
