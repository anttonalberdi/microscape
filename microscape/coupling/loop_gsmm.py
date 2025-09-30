from __future__ import annotations
from pathlib import Path
from typing import Dict, Tuple, Callable
import yaml
import numpy as np

from ..utils.resources import get_packaged_path

ProgressCB = Callable[[int, int], None]

def load_config(yaml_path: str) -> dict:
    """Load YAML config; if relative/not found, try packaged examples."""
    p = Path(yaml_path)
    if p.exists():
        return yaml.safe_load(p.read_text())

    # try inside package (e.g., microscape/examples/demo2/community_sbml.yml)
    packaged = get_packaged_path(yaml_path)
    pp = Path(packaged)
    if pp.exists():
        return yaml.safe_load(pp.read_text())

    raise FileNotFoundError(f"Config not found: {yaml_path}")

def _resolve_relative_to_config(path_like: str, config_dir: Path) -> Path:
    """Resolve a model path robustly: absolute -> as is; else relative to config dir; else packaged."""
    p = Path(path_like)
    if p.is_absolute() and p.exists():
        return p

    # relative to the YAML’s directory
    cand = (config_dir / p).resolve()
    if cand.exists():
        return cand

    # packaged fallback
    packaged = get_packaged_path(str(p))
    pp = Path(packaged)
    if pp.exists():
        return pp

    raise FileNotFoundError(f"SBML model not found: {path_like} (looked in {config_dir})")

# ---- Your existing demo core (minimal placeholder below) ----

def run_demo_gsmm(config_path: str, outdir: str, progress_cb: ProgressCB | None = None) -> Tuple[Dict[str, np.ndarray], Dict[str, float]]:
    """
    Minimal scaffold: resolves model paths and then runs your existing
    ABM+dFBA+RD loop. Replace the 'fake sim' with your actual code.
    """
    cfg = load_config(config_path)
    config_dir = Path(get_packaged_path(config_path)).parent if not Path(config_path).exists() else Path(config_path).resolve().parent

    # Normalize model paths (support both dict and list styles if you ever change schema)
    models_cfg = cfg.get("models", {})
    for k, raw in list(models_cfg.items()):
        models_cfg[k] = str(_resolve_relative_to_config(raw, config_dir))

    # ---- TODO: call your actual simulation using 'models_cfg' paths ----
    H = W = int(cfg.get("grid", {}).get("size_px", 128))
    steps = int(cfg.get("grid", {}).get("steps", 100))
    if progress_cb:
        progress_cb(1, 100)

    # FAKE fields just to keep the interface (replace with real outputs)
    fields = {
        "glc": np.zeros((H, W), dtype=float),
        "lac": np.zeros((H, W), dtype=float),
        "ac":  np.zeros((H, W), dtype=float),
        "but": np.zeros((H, W), dtype=float),
    }
    # Simple “progress” ticks
    for i in range(2, 101, max(1, 100 // max(1, steps))):
        if progress_cb:
            progress_cb(i, 100)

    summary = {
        "n_steps": float(steps),
        "grid_size_px": float(H),
        "models_resolved": float(len(models_cfg)),
    }
    return fields, summary