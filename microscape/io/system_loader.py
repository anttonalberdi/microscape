# microscape/io/system_loader.py
from __future__ import annotations
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import yaml

# ---------- helpers ----------

def _read_yaml(p: Path) -> dict:
    with p.open("r", encoding="utf-8") as fh:
        return yaml.safe_load(fh) or {}

def _resolve(base: Path, maybe: Optional[str]) -> Optional[Path]:
    if not maybe:
        return None
    p = Path(maybe)
    return p if p.is_absolute() else (base / p)

# ---------- public API ----------

def load_system(system_yml: Path) -> dict:
    """
    Load system.yml and resolve:
      - root
      - paths (config_dir, environments_dir, spots_dir, microbes_dir)
      - ecology_cfg, metabolism_cfg (resolved to files)
      - environment_files (list of Paths)
      - microbe_registry (dict id -> microbe.yml Path)
    """
    system_yml = Path(system_yml).resolve()
    root = system_yml.parent
    sysd = _read_yaml(system_yml).get("system", {})

    # paths section (with defaults)
    paths = sysd.get("paths") or {}
    config_dir = _resolve(root, paths.get("config_dir")) or (root / "config")
    envs_dir   = _resolve(root, paths.get("environments_dir")) or (root / "environments")
    spots_dir  = _resolve(root, paths.get("spots_dir")) or (root / "spots")
    microbes_dir = _resolve(root, paths.get("microbes_dir")) or (root / "microbes")

    # config.ecology / config.metabolism (relative to config_dir unless absolute)
    config_map = sysd.get("config") or {}
    ecology_cfg    = _resolve(config_dir, config_map.get("ecology")) if config_map.get("ecology") else None
    metabolism_cfg = _resolve(config_dir, config_map.get("metabolism")) if config_map.get("metabolism") else None

    # environments registry: supports [{id,file}, ...] or ["E001.yml", "E002.yml"] or ["E001", ...]
    env_specs = ((sysd.get("registry") or {}).get("environments") or [])
    env_files: List[Path] = []
    if env_specs:
        for it in env_specs:
            if isinstance(it, str):
                f = it if it.endswith(".yml") else f"{it}.yml"
                env_files.append((envs_dir / f).resolve())
            elif isinstance(it, dict):
                fid = it.get("file") or (it.get("id") and f"{it['id']}.yml")
                if fid:
                    p = Path(fid)
                    env_files.append(p if p.is_absolute() else (envs_dir / p).resolve())
    else:
        env_files = sorted(envs_dir.glob("*.yml"))

    env_files = [p for p in env_files if p.exists()]

    # microbes registry: [{id,file}, ...] or ["M0001.yml"] or ["M0001", ...]
    micro_specs = ((sysd.get("registry") or {}).get("microbes") or [])
    microbe_registry: Dict[str, Path] = {}
    if micro_specs:
        for it in micro_specs:
            if isinstance(it, str):
                mid = it.removesuffix(".yml")
                f = (microbes_dir / f"{mid}.yml").resolve()
                microbe_registry[mid] = f
            elif isinstance(it, dict):
                mid = it.get("id")
                mf  = it.get("file") or (mid and f"{mid}.yml")
                if mid and mf:
                    p = Path(mf)
                    microbe_registry[mid] = (p if p.is_absolute() else (microbes_dir / p)).resolve()
    else:
        # auto-discover
        for p in sorted(microbes_dir.glob("*.yml")):
            microbe_registry[p.stem] = p.resolve()

    return {
        "root": str(root),
        "system": sysd,
        "paths": {
            "config_dir": str(config_dir),
            "environments_dir": str(envs_dir),
            "spots_dir": str(spots_dir),
            "microbes_dir": str(microbes_dir),
        },
        "ecology_cfg": str(ecology_cfg) if ecology_cfg else None,
        "metabolism_cfg": str(metabolism_cfg) if metabolism_cfg else None,
        "environment_files": env_files,
        "microbe_registry": microbe_registry,
    }

def iter_spot_files_for_env(env_file: Path, sys_paths: Dict) -> List[Tuple[str, Path]]:
    """Return [(spot_id, spot_path)] for an environment.yml, honoring per-env list or dir."""
    env = _read_yaml(env_file).get("environment", {})
    base = env_file.parent

    # Resolve candidate spots base directory:
    env_spots_dir = env.get("spots_dir")
    if env_spots_dir:
        spots_base = _resolve(base, env_spots_dir)
    else:
        # fall back to system-level spots_dir colocated with system root
        sys_spots_dir = sys_paths.get("spots_dir")
        spots_base = Path(sys_spots_dir) if sys_spots_dir else (base / "spots")

    out: List[Tuple[str, Path]] = []

    # explicit list
    listed = env.get("spots")
    if isinstance(listed, list) and listed:
        for s in listed:
            sid = s.get("id") or s.get("name")
            f = s.get("file")
            if not sid or not f:
                continue
            p = Path(f)
            spath = p if p.is_absolute() else ((spots_base / p) if p.parent == Path(".") else (base / p))
            out.append((sid, spath.resolve()))
        return out

    # else glob directory
    if spots_base and spots_base.exists():
        for p in sorted(spots_base.glob("*.yml")):
            out.append((p.stem, p.resolve()))
    return out

def read_spot_yaml(spot_path: Path) -> dict:
    return _read_yaml(Path(spot_path))

def read_microbe_yaml(microbe_id_or_path: str | Path, sys_info: dict | None = None) -> Optional[dict]:
    """
    Load microbe YAML either by absolute/relative path, or by ID via sys_info['microbe_registry'].
    """
    if isinstance(microbe_id_or_path, (str, Path)):
        p = Path(microbe_id_or_path)
        if p.suffix.lower() == ".yml" and p.exists():
            return _read_yaml(p)
        # try resolve via registry
        if sys_info and isinstance(microbe_id_or_path, str):
            reg = sys_info.get("microbe_registry", {})
            mp = reg.get(microbe_id_or_path)
            if mp and Path(mp).exists():
                return _read_yaml(Path(mp))
    return None
