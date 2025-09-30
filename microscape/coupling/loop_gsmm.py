from __future__ import annotations
from pathlib import Path
from typing import Dict, Tuple
import numpy as np
import yaml
import importlib.resources as ir


from ..rd.simple import diffuse_step
from ..io.sbml import load_sbml_model, map_exchanges, require_cobra
from ..metabolism.fba import run_fba_with_bounds

def load_config(yaml_path: str) -> dict:
    p = Path(yaml_path)
    if not p.exists():
        # Try to load from packaged resources
        # First assume path is relative to packaged examples
        candidates = []
        # packaged examples under top-level "examples"
        candidates.append(("examples/00_synthetic/community_sbml.yml", None))
        # if user passed a relative path, try it under packaged root
        if not yaml_path.startswith("examples/"):
            candidates.append((yaml_path, None))

        for rel, _ in candidates:
            try:
                with ir.files("microscape").joinpath(rel).open("rb") as fh:
                    import yaml as _yaml
                    return _yaml.safe_load(fh.read())
            except Exception:
                pass
        # As a last attempt, look next to this module (installed tree)
        pkg_root = Path(__file__).resolve().parents[2]  # .../microscape
        alt = pkg_root / "examples" / "00_synthetic" / "community_sbml.yml"
        if alt.exists():
            import yaml as _yaml
            return _yaml.safe_load(alt.read_text())

        raise FileNotFoundError(f"Config not found: {p} (and no packaged fallback found)")
    import yaml
    return yaml.safe_load(p.read_text())

def init_fields(cfg: dict, H: int, W: int) -> Dict[str, np.ndarray]:
    fields = {}
    init = cfg.get("fields", {}).get("init", {})
    for name, val in init.items():
        arr = np.full((1, H, W), float(val), dtype=float)
        fields[name] = arr
    # synthetic fibre patch converted to glucose via fibre_to_seeds.glc_factor
    env = cfg.get("environment", {})
    patch = env.get("fibre_patch", {"shape": "square", "size_px": 10})
    size = int(patch.get("size_px", 10))
    y0, x0 = H//2 - size//2, W//2 - size//2
    mask = np.zeros((1, H, W), dtype=float)
    mask[:, y0:y0+size, x0:x0+size] = 1.0
    fibre_to = env.get("fibre_to_seeds", {})
    glc_fac = float(fibre_to.get("glc_factor", 1.0))
    if "glc" in fields:
        fields["glc"] = fields["glc"] + glc_fac * mask
    else:
        fields["glc"] = glc_fac * mask
    return fields

def dC_from_flux(flux_mmol_g_h: float, alpha: float, dt_s: float) -> float:
    return alpha * flux_mmol_g_h * dt_s

def run_demo_gsmm(config_path: str, outdir: str, progress_cb=None) -> Tuple[Dict[str, np.ndarray], Dict]:
    cfg = load_config(config_path)
    H = W = int(cfg.get("grid", {}).get("size_px", 128))
    dt = float(cfg.get("grid", {}).get("dt_s", 2.0))
    steps = int(cfg.get("grid", {}).get("steps", 200))
    voxel_um = float(cfg.get("grid", {}).get("voxel_um", 10))

    D = cfg.get("fields", {}).get("diffusion_um2_s", {})
    decay = cfg.get("fields", {}).get("decay_s", {})
    alpha = float(cfg.get("alpha_flux_to_dC", 1e-3))

    fields = init_fields(cfg, H, W)
    for k in D:
        fields.setdefault(k, np.zeros((1, H, W), dtype=float))

    mucosa_px = int(cfg.get("environment", {}).get("mucosa_band_px", 5))
    mucosa_mask = np.zeros((1, H, W), dtype=bool)
    mucosa_mask[:, 0:mucosa_px, :] = True

    require_cobra()
    model_paths = cfg.get("models", {})
    if not model_paths:
        raise RuntimeError("No models defined in config under 'models'.")
    models = {name: load_sbml_model(path) for name, path in model_paths.items()}

    exmap_cfg = cfg.get("fields_to_exchanges") or {}
    if not exmap_cfg:
        exmap_file = Path(__file__).resolve().parent.parent / "config" / "exchanges.yml"
        if exmap_file.exists():
            import yaml as _y; exmap_cfg = _y.safe_load(exmap_file.read_text()).get("fields_to_exchanges", {})
    if not exmap_cfg:
        raise RuntimeError("No exchange mapping provided (fields_to_exchanges)." )

    # main loop
    for i in range(steps):
        mean_conc = {k: float(np.nanmean(v)) for k, v in fields.items()}

        kmap = cfg.get("bounds", {}).get("uptake_k", {})
        umax = cfg.get("bounds", {}).get("max_ub", {})
        bounds = {}
        for fname, rxn_id in exmap_cfg.items():
            kfac = float(kmap.get(fname, 0.0))
            ub = min(float(umax.get(fname, 10.0)), kfac * mean_conc.get(fname, 0.0))
            bounds[rxn_id] = ub

        sec = {f: 0.0 for f in fields.keys()}
        for mname, path in model_paths.items():
            model = models[mname]
            sol = run_fba_with_bounds(model, bounds, use_pfba=True)
            flx = sol.fluxes
            for fname, rxn_id in exmap_cfg.items():
                if rxn_id in flx.index:
                    v = float(flx[rxn_id])
                    sec[fname] = sec.get(fname, 0.0) + v

        for fname, net_flux in sec.items():
            if fname not in fields:
                continue
            dC = dC_from_flux(net_flux, alpha=alpha, dt_s=dt)
            fields[fname] = fields[fname] + dC

        for fname, arr in fields.items():
            Di = float(D.get(fname, 0.0))
            dec = float(decay.get(fname, 0.0))
            if Di > 0 or dec > 0:
                fields[fname] = diffuse_step(arr, D=Di, dt=dt, dx=voxel_um, decay=dec)

        if progress_cb is not None:
            progress_cb(i+1, steps)

    summary = {}
    for k, v in fields.items():
        a2d = v[0]
        summary[f"mean_{k}"] = float(np.nanmean(a2d))
        summary[f"max_{k}"] = float(np.nanmax(a2d))
        if k in ("but","butyrate","but_e") :
            summary["scfa_at_mucosa"] = float(np.nanmean(a2d[mucosa_mask[0]]))

    fields["mucosa_mask"] = mucosa_mask.astype(float)
    return fields, summary
