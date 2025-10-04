from __future__ import annotations
from pathlib import Path
from typing import Any, Dict, Callable, List
import yaml, json, csv

def _read_yaml(p: Path) -> dict:
    return yaml.safe_load(p.read_text())

def _write_yaml(p: Path, obj: dict):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(yaml.safe_dump(obj, sort_keys=False))

# --- scoring helpers ---
def _score_clip_linear(expr: float, thr: float) -> float:
    if thr <= 0: return 0.0
    s = expr / thr
    return 0.0 if s < 0 else (1.0 if s > 1 else s)

def _score_binary(expr: float, thr: float) -> float:
    return 1.0 if expr >= thr else 0.0

def _score_logistic(expr: float, thr: float, k: float = 1.0) -> float:
    import math
    return 1.0 / (1.0 + math.exp(-k * (expr - thr)))

def _get_scoring_fn(cfg: dict) -> Callable[[float, float], float]:
    mode = ((cfg.get("ecology") or {}).get("scoring") or {}).get("mode", "clip_linear")
    if mode == "clip_linear":
        return _score_clip_linear
    if mode == "binary":
        return _score_binary
    if mode == "logistic":
        pars = ((cfg.get("ecology") or {}).get("scoring") or {}).get("logistic", {})
        k = float(pars.get("k", 1.0))
        return lambda expr, thr: _score_logistic(expr, thr, k)
    raise ValueError(f"Unknown scoring mode: {mode}")

def _infer_traits_for_microbe(m_id: str, gene_expr: Dict[str, float], traits: List[dict], cfg_scoring: dict) -> Dict[str, Any]:
    out = {}
    scoring_fn = _get_scoring_fn({"ecology": cfg_scoring})
    for tr in traits:
        tr_id   = tr["id"]
        states  = tr["states"]
        pos     = tr.get("positive_state", states[0])
        # overrides
        ov      = (tr.get("overrides") or {}).get(m_id, {})
        markers = ov.get("markers", tr.get("markers", []))
        thr     = float(ov.get("threshold", tr.get("threshold", 0.0)))
        # max over marker TPMs
        expr = max(float(gene_expr.get(g, 0.0)) for g in markers) if markers else 0.0
        score = scoring_fn(expr, thr) if thr > 0 else 0.0
        # assign state (two-state; first is positive unless overridden)
        neg = states[1] if pos == states[0] and len(states) > 1 else (states[0] if len(states) > 1 else states[0])
        state = pos if score >= 0.5 else neg
        out[tr_id] = {
            "state": state,
            "score": round(float(score), 3),
            "marker_max_TPM": round(float(expr), 3),
            "markers": markers,
            "threshold": thr,
        }
    return out

def _enrich_spot_yaml(spot_yaml: Path, ecology_cfg: dict) -> dict:
    spot = _read_yaml(spot_yaml)
    meas = (spot.get("measurements") or {})
    tblock = (meas.get("transcripts") or {})
    schema = tblock.get("schema", "genes")
    if schema != "genes":
        raise ValueError(f"{spot_yaml.name}: transcripts.schema must be 'genes', got '{schema}'")
    per_microbe = tblock.get("values", {})
    traits = ((ecology_cfg.get("ecology") or {}).get("traits") or [])
    scoring = ((ecology_cfg.get("ecology") or {}).get("scoring") or {})
    # build ecology states per microbe
    eco = {}
    for m_id, gene_expr in per_microbe.items():
        eco[m_id] = _infer_traits_for_microbe(m_id, gene_expr, traits, scoring)
    # write into structure
    meas["ecology"] = {
        "schema": "traits",
        "traits": [t["id"] for t in traits],
        "states": eco
    }
    spot["measurements"] = meas
    return spot

def _iter_environments(system_path: Path) -> List[Path]:
    sysd = _read_yaml(system_path)
    reg = ((sysd.get("system") or {}).get("registry") or {})
    env_ids = reg.get("environments") or []
    # default: environments/ENV_ID.yml next to system.yml
    base = system_path.parent / "environments"
    env_files = []
    for eid in env_ids:
        p = base / f"{eid}.yml"
        if p.exists():
            env_files.append(p)
    return env_files

def _iter_spot_files(env_yaml: Path) -> List[Path]:
    envd = _read_yaml(env_yaml)
    spots = ((envd.get("environment") or {}).get("spots") or [])
    res = []
    for s in spots:
        sf = s.get("file")
        if not sf: continue
        p = (env_yaml.parent / sf).resolve()
        if p.exists():
            res.append(p)
    return res

def profile_system(system_yml: Path, ecology_rules_yml: Path, outdir: Path,
                   overwrite: bool = False, progress_cb=None) -> Dict[str, Any]:
    ecology_cfg = _read_yaml(ecology_rules_yml)
    env_files = _iter_environments(system_yml)

    n_env = len(env_files)
    n_spots = 0
    out_spots = 0

    # mirror structure under outdir: environments/ENV/spots/SPOT.yml
    for ei, env_path in enumerate(env_files):
        envd = _read_yaml(env_path)
        eid = (envd.get("environment") or {}).get("id", env_path.stem)
        spot_files = _iter_spot_files(env_path)
        n_spots += len(spot_files)
        for si, sp in enumerate(spot_files):
            enriched = _enrich_spot_yaml(sp, ecology_cfg)
            rel = sp.relative_to(env_path.parent)
            out_p = outdir / "environments" / eid / rel
            if out_p.exists() and not overwrite:
                # keep existing; but still count
                pass
            else:
                _write_yaml(out_p, enriched)
            out_spots += 1
            if progress_cb:
                # rough progress: per spot
                done = (out_spots / max(1, n_spots)) * 100.0
                progress_cb(done)

    # Build a tiny summary table (CSV & JSON)
    # For speed, re-read enriched spots under outdir
    rows = []
    for env_dir in (outdir / "environments").glob("*"):
        env_yaml = env_dir / f"{env_dir.name}.yml"  # may not exist in outdir; skip
        for sp in env_dir.glob("spots/*.yml"):
            s = _read_yaml(sp)
            env_id = env_dir.name
            spot_id = (s.get("spot") or {}).get("name", sp.stem)
            eco = ((s.get("measurements") or {}).get("ecology") or {})
            states = eco.get("states", {})
            # flatten: for each microbe and trait, store state & score
            for m_id, trmap in states.items():
                for tr_id, rec in trmap.items():
                    rows.append({
                        "environment": env_id,
                        "spot": spot_id,
                        "microbe": m_id,
                        "trait": tr_id,
                        "state": rec.get("state"),
                        "score": rec.get("score"),
                    })

    # CSV text
    from io import StringIO
    sio = StringIO()
    w = csv.DictWriter(sio, fieldnames=["environment","spot","microbe","trait","state","score"])
    w.writeheader()
    for r in rows: w.writerow(r)
    summary_csv = sio.getvalue()

    summary_json = {
        "n_environments": n_env,
        "n_spots_processed": out_spots,
        "n_rows": len(rows),
        "traits": [t["id"] for t in ((ecology_cfg.get("ecology") or {}).get("traits") or [])]
    }
    return {"summary_csv": summary_csv, "summary_json": summary_json}
