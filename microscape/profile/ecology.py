from __future__ import annotations
from typing import Dict, Any, List
from pathlib import Path
import yaml

def load_rules(p: Path) -> Dict[str, Any]:
    return yaml.safe_load(Path(p).read_text())

def _spot_code(sd: dict, spot_yml: Path) -> str:
    # Prefer 'name', then 'id'; fallback to filename stem
    return (sd.get("name")
            or sd.get("id")
            or spot_yml.stem)

def score_clip_linear(x: float, thr: float, max_val: float | None = None) -> float:
    # Simple saturation: 0 below thr/2, ramp to 1 at thr, then clip
    if x is None or x <= 0:
        return 0.0
    lo = 0.5 * thr
    hi = thr
    if x <= lo:
        return 0.0
    if x >= hi:
        return 1.0
    return (x - lo) / (hi - lo)

def profile_spot(spot_yml: Path, rules: Dict[str, Any]) -> List[Dict[str, Any]]:
    sd = yaml.safe_load(spot_yml.read_text())["spot"]

    # robust spot code
    spot_code = _spot_code(sd, spot_yml)

    # microbes block
    microbes = (sd.get("measurements", {}).get("microbes") or {})
    mvals = microbes.get("values") or {}  # {microbe_id: abundance}

    # transcripts block (can be missing or empty)
    transcripts = (sd.get("measurements", {}).get("transcripts") or {})
    unit = transcripts.get("unit") or "expr"
    # We support either `per_microbe:` or `values:` (dict of dicts)
    per_microbe = transcripts.get("per_microbe") or transcripts.get("values") or {}

    traits = (rules.get("ecology", {}).get("traits") or [])
    scoring_mode = (rules.get("ecology", {}).get("scoring", {}) or {}).get("mode", "clip_linear")

    rows: List[Dict[str, Any]] = []
    for m_id, abundance in mvals.items():
        gexp: Dict[str, float] = per_microbe.get(m_id, {}) if isinstance(per_microbe, dict) else {}
        out: Dict[str, Any] = {
            "spot": spot_code,
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
            mean_expr = sum(expr_vals) / len(expr_vals) if expr_vals else 0.0

            if scoring_mode == "clip_linear":
                score = score_clip_linear(mean_expr, thr)
            else:
                score = 1.0 if mean_expr >= thr else 0.0

            states = t.get("states", ["on", "off"])
            on_state = states[0] if states else "on"
            off_state = states[1] if len(states) > 1 else "off"

            out[f"{tid}_state"] = on_state if score >= 0.5 else off_state
            out[f"{tid}_score"] = round(score, 3)
            out[f"{tid}_expr_{unit}"] = round(mean_expr, 3)

        rows.append(out)

    return rows
