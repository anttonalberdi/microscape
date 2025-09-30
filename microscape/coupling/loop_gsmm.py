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

def _seed_fibre_patch(fields, cfg, H, W):
    patch = (cfg.get("environment", {}) or {}).get("fibre_patch", {})
    seeds = (cfg.get("environment", {}) or {}).get("fibre_to_seeds", {})
    if not patch or "glc_value" not in seeds:
        return
    size = int(patch.get("size_px", 12))
    cx_rel, cy_rel = patch.get("center", [0.5, 0.5])
    cx, cy = int(cx_rel * W), int(cy_rel * H)
    x0, x1 = max(0, cx - size//2), min(W, cx + size//2)
    y0, y1 = max(0, cy - size//2), min(H, cy + size//2)
    fields["glc"][y0:y1, x0:x1] = float(seeds["glc_value"])


def run_demo_gsmm(config_path: str, outdir: str, progress_cb: ProgressCB | None = None) -> Tuple[Dict[str, np.ndarray], Dict[str, float]]:
    """
    Demo 2: toy cross-feeding simulation (FD, LU, BP) with reaction–diffusion.
    No GSMM solver yet; uses simple local kinetics + diffusion to produce visible patterns.
    """
    cfg = load_config(config_path)
    # Resolve model paths so the demo validates that files exist (forward-compatible)
    config_dir = Path(get_packaged_path(config_path)).parent if not Path(config_path).exists() else Path(config_path).resolve().parent
    for k, raw in list((cfg.get("models", {}) or {}).items()):
        _ = _resolve_relative_to_config(raw, config_dir)  # just to check; not parsed here

    # Grid/time
    H = W = int(cfg.get("grid", {}).get("size_px", 128))
    dt = float(cfg.get("grid", {}).get("dt_s", 2.0))
    steps = int(cfg.get("grid", {}).get("steps", 200))
    voxel_um = float(cfg.get("grid", {}).get("voxel_um", 10.0))
    h = voxel_um  # µm per pixel

    # Diffusion (µm^2/s)
    diffs = cfg.get("fields", {}).get("diffusion_um2_s", {}) or {}
    D_glc = float(diffs.get("glc", 600))
    D_lac = float(diffs.get("lac", 500))
    D_ac  = float(diffs.get("ac",  500))
    D_but = float(diffs.get("but", 400))

    # Decay (1/s)
    decays = cfg.get("fields", {}).get("decay_s", {}) or {}
    k_decay_glc = float(decays.get("glc", 0.0))
    k_decay_lac = float(decays.get("lac", 0.0))
    k_decay_ac  = float(decays.get("ac",  0.0))
    k_decay_but = float(decays.get("but", 0.0))

    # Flux → ∆C scale (a.u. per mmol/gDW/h; just a demo scale factor)
    alpha = float(cfg.get("alpha_flux_to_dC", 5e-3))

    # Reaction rates (1/s per concentration unit) – use YAML bounds as “uptake_k” if present
    uptake_k = (cfg.get("bounds", {}) or {}).get("uptake_k", {}) or {}
    k_glc = float(uptake_k.get("glc", 0.05))
    k_lac = float(uptake_k.get("lac", 0.04))
    k_ac  = float(uptake_k.get("ac",  0.03))
    # Yields
    y_fd_lac, y_fd_ac  = 0.6, 0.4
    y_lu_but           = 0.9
    y_bp_but           = 0.8

    # Fields
    glc = np.zeros((H, W), dtype=float)
    lac = np.zeros((H, W), dtype=float)
    ac  = np.zeros((H, W), dtype=float)
    but = np.zeros((H, W), dtype=float)

    fields = {"glc": glc, "lac": lac, "ac": ac, "but": but}

    # Seed fibre glucose patch from YAML
    _seed_fibre_patch(fields, cfg, H, W)

    def laplacian(A: np.ndarray, D: float) -> np.ndarray:
        # 5-point laplacian with no-flux (Neumann) boundaries via edge reflection
        # Pad by 1 with edge values
        A_p = np.pad(A, 1, mode="edge")
        lap = (A_p[1:-1,2:] + A_p[1:-1,:-2] + A_p[2:,1:-1] + A_p[:-2,1:-1] - 4*A) / (h*h)
        return D * lap

    # crude stability check (not strict, just a guard)
    for D in (D_glc, D_lac, D_ac, D_but):
        if dt * D / (h*h) > 0.25:
            # reduce dt if wildly unstable
            dt = 0.25 * (h*h) / max(D, 1e-9)

    # Time loop
    for t in range(steps):
        # Local “kinetics” (mass-action like). Use alpha to push visible signal.
        r_fd = alpha * k_glc * glc            # glucose consumption rate
        r_lu = alpha * k_lac * lac            # lactate consumption
        r_bp = alpha * k_ac  * ac             # acetate consumption

        # Apply reactions
        glc_new = glc - r_fd * dt
        lac_new = lac + (y_fd_lac * r_fd - r_lu) * dt
        ac_new  = ac  + (y_fd_ac  * r_fd - r_bp) * dt
        but_new = but + (y_lu_but * r_lu + y_bp_but * r_bp) * dt

        # Decay
        if k_decay_glc: glc_new -= k_decay_glc * glc * dt
        if k_decay_lac: lac_new -= k_decay_lac * lac * dt
        if k_decay_ac:  ac_new  -= k_decay_ac  * ac  * dt
        if k_decay_but: but_new -= k_decay_but * but * dt

        # Clamp non-negativity
        glc_new = np.clip(glc_new, 0, None)
        lac_new = np.clip(lac_new, 0, None)
        ac_new  = np.clip(ac_new,  0, None)
        but_new = np.clip(but_new, 0, None)

        # Diffusion step
        glc = glc_new + dt * laplacian(glc_new, D_glc)
        lac = lac_new + dt * laplacian(lac_new, D_lac)
        ac  = ac_new  + dt * laplacian(ac_new,  D_ac)
        but = but_new + dt * laplacian(but_new, D_but)

        # prepare for next
        fields["glc"], fields["lac"], fields["ac"], fields["but"] = glc, lac, ac, but

        if progress_cb and (t % max(1, steps // 100) == 0):
            progress_cb(int(100 * (t + 1) / steps), 100)

    # Simple summary
    summary = {
        "n_steps": float(steps),
        "grid_size_px": float(H),
        "mean_glc": float(glc.mean()),
        "mean_lac": float(lac.mean()),
        "mean_ac":  float(ac.mean()),
        "mean_but": float(but.mean()),
    }
    return fields, summary
