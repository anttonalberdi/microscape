from __future__ import annotations
from typing import Dict, Any, List
from pathlib import Path
import yaml, csv, json

def load_rules(p: Path) -> Dict[str, Any]:
    return yaml.safe_load(Path(p).read_text())

def score_clip_linear(x: float, thr: float, max_val: float = None) -> float:
    # simple saturation: 0 below thr/2, ramp to 1 at thr, then clip
    if x is None:
        return 0.0
    if x <= 0:
        return 0.0
    # ramp around threshold for a bit of robustness
    lo = 0.5 * thr
    hi = thr
    if x <= lo:
        return 0.0
    if x >= hi:
        return 1.0
    return (x - lo) / (hi - lo)

def profile_spot(spot_yml: Path, rules: Dict[str, Any]) -> List[Dict[str, Any]]:
    sd = yaml.safe_load(spot_yml.read_text())["spot"]
    microbes = (sd.get("measurements", {}).get("microbes", {}) or {})
    mvals = microbes.get("values") or {}
    transcripts = (sd.get("measurements", {}).get("transcripts") or {})
    unit = transcripts.get("unit")
    schema = transcripts.get("schema")
    per = transcripts.get("per_microbe") or transcripts.get("values") or {}
    traits = (rules.get("ecology", {}).get("traits") or [])
    scoring_mode = (rules.get("ecology", {}).get("scoring", {}) or {}).get("mode","clip_linear")

    rows = []
    for m_id, abundance in mvals.items():
        gexp: Dict[str, float] = per.get(m_id, {}) if isinstance(per, dict) else {}
        out: Dict[str, Any] = {
            "spot": sd.get("name"),
            "microbe": m_id,
            "abundance": abundance,
        }
        # compute trait states & scores
        for t in traits:
            tid = t["id"]
            markers = t.get("markers", [])
            thr = float(t.get("threshold", 10))
            # combine marker expression (e.g., mean of listed genes)
            expr_vals = [float(gexp.get(g, 0.0)) for g in markers]
            mean_expr = sum(expr_vals)/len(expr_vals) if expr_vals else 0.0
            if scoring_mode == "clip_linear":
                score = score_clip_linear(mean_expr, thr)
            else:
                score = 1.0 if mean_expr >= thr else 0.0
            states = t.get("states", ["on","off"])
            # pick “on” state if score >= 0.5
            out[f"{tid}_state"] = states[0] if score >= 0.5 else (states[1] if len(states)>1 else "off")
            out[f"{tid}_score"] = round(score, 3)
            out[f"{tid}_expr_{unit or 'expr'}"] = round(mean_expr, 3)
        rows.append(out)
    return rows
