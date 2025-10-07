# microscape/io/system_loader.py
from __future__ import annotations
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import yaml

def _read_yaml(p: Path) -> dict:
    with p.open("r", encoding="utf-8") as fh:
        return yaml.safe_load(fh) or {}

def _resolve(base: Path, maybe: Optional[str]) -> Optional[Path]:
    if not maybe: return None
    p = Path(maybe)
    return (base / p) if not p.is_absolute() else p

def load_system(system_yml: Path) -> dict:
    system_yml = Path(system_yml).resolve()
    root = system_yml.parent
    sysd = _read_yaml(system_yml)["system"]

    paths: Dict = sysd.get("paths") or {}
    envs_dir = _resolve(root, paths.get("environments_dir")) or (root / "environments")
    config_dir = _resolve(root, paths.get("config_dir")) or (root / "config")
    spots_dir  = _resolve(root, paths.get("spots_dir")) or (root / "spots")
    microbes_dir = _resolve(root, paths.get("microbes_dir")) or (root / "microbes")

    ecology_rel = (sysd.get("config") or {}).get("ecology")
    metab_rel   = (sysd.get("config") or {}).get("metabolism")
    ecology_cfg = _resolve(config_dir, ecology_rel) if ecology_rel else None
    metab_cfg   = _resolve(config_dir, metab_rel) if metab_rel else None

    env_specs = ((sysd.get("registry") or {}).get("environments") or [])
    microbe_specs = ((sysd.get("registry") or {}).get("microbes") or [])

    env_files: List[Path] = []
    for item in env_specs:
        if isinstance(item, str):
            env_files.append(_resolve(envs_dir, item if item.endswith(".yml") else f"{item}.yml"))
        elif isinstance(item, dict):
            fid = item.get("file") or f"{item.get('id')}.yml"
            env_files.append(_resolve(envs_dir, fid))
    env_files = [p for p in env_files if p and p.exists()]

    microbe_files: Dict[str, Path] = {}
    for item in microbe_specs:
        if not isinstance(item, dict): continue
        mid = item.get("id"); mf = item.get("file")
        if not mid or not mf: continue
        mp = _resolve(microbes_dir, mf)
        if mp and mp.exists():
            microbe_files[mid] = mp

    return {
        "root": root,
        "system": sysd,
        "paths": {
            "envs_dir": envs_dir, "config_dir": config_dir,
            "spots_dir": spots_dir, "microbes_dir": microbes_dir
        },
        "configs": {"ecology": ecology_cfg, "metabolism": metab_cfg},
        "environment_files": env_files,
        "microbe_files": microbe_files,
    }

def iter_spot_files_for_env(env_file: Path, sys_paths: Dict) -> List[Tuple[str, Path]]:
    env = _read_yaml(env_file)["environment"]
    base = env_file.parent
    env_spots_dir = env.get("spots_dir")
    sys_spots_dir = sys_paths.get("spots_dir")
    if env_spots_dir:
        spots_base = _resolve(base, env_spots_dir)
    elif sys_spots_dir:
        spots_base = sys_spots_dir
    else:
        spots_base = base / "spots"

    out: List[Tuple[str, Path]] = []
    if env.get("spots"):
        for s in env["spots"]:
            sid = s.get("id"); f = s.get("file")
            if not sid or not f: continue
            p = Path(f)
            spath = (spots_base / p) if (p.parent == Path(".")) else _resolve(base, f)
            out.append((sid, spath))
        return out

    if spots_base and spots_base.exists():
        for p in sorted(spots_base.glob("*.yml")):
            out.append((p.stem, p))
    return out

def read_spot_yaml(spot_path: Path) -> dict:
    return _read_yaml(spot_path).get("spot", {})

def read_microbe_yaml(mid: str, sys_info: dict) -> Optional[dict]:
    p = sys_info.get("microbe_files", {}).get(mid)
    if not p: return None
    d = _read_yaml(p).get("microbe", {})
    if not d: return None
    d["_file"] = p
    return d

def spot_transcripts(spot: dict) -> Dict[str, Dict[str, float]]:
    """Return {microbe_id: {gene_id: TPM}} from a spot dict."""
    tr = ((spot.get("measurements") or {}).get("transcripts") or {})
    vals = tr.get("values") or {}
    out = {}
    for mid, gmap in vals.items():
        if not isinstance(gmap, dict): continue
        out[str(mid)] = {str(g): float(v) for g, v in gmap.items()}
    return out

def spot_metabolites(spot: dict) -> Dict[str, float]:
    mets = ((spot.get("measurements") or {}).get("metabolites") or {})
    vals = mets.get("values") or {}
    return {str(k): float(v) for k, v in vals.items()}
