# microscape/io/system_loader.py
from __future__ import annotations
from pathlib import Path
from typing import List, Dict, Tuple, Any, Optional
import yaml

__all__ = [
    "load_system",
    "iter_spot_files_for_env",
    "read_spot_yaml",
    "read_microbe_yaml",
    "load_microbe_registry",
]

def _read_yaml(p: Path) -> dict:
    with p.open("r", encoding="utf-8") as fh:
        return yaml.safe_load(fh) or {}

def _resolve(base: Path, maybe: Optional[str]) -> Optional[Path]:
    if not maybe:
        return None
    p = Path(maybe)
    return (base / p) if not p.is_absolute() else p

def load_system(system_yml: Path) -> dict:
    """
    Load system.yml and resolve key directories/files relative to the system root.
    Expected structure:
      system.paths.{config_dir,environments_dir,spots_dir,microbes_dir}
      system.config.ecology / metabolism
      system.registry.environments: [{id, file} | "E001" | "E001.yml"]
      system.registry.microbes:     [{id, file} | "M0001" | "M0001.yml"]
    """
    system_yml = Path(system_yml).resolve()
    root = system_yml.parent
    sysd = _read_yaml(system_yml).get("system") or {}

    paths: Dict[str, Any] = sysd.get("paths") or {}
    envs_dir = _resolve(root, paths.get("environments_dir")) or (root / "environments")
    config_dir = _resolve(root, paths.get("config_dir")) or (root / "config")
    # Keep for callers that need it
    paths.setdefault("environments_dir", str(envs_dir))
    paths.setdefault("config_dir", str(config_dir))
    if "spots_dir" in paths:
        paths["spots_dir"] = str(_resolve(root, paths["spots_dir"]))
    if "microbes_dir" in paths:
        paths["microbes_dir"] = str(_resolve(root, paths["microbes_dir"]))

    # Resolve config files (ecology/metabolism)
    cfg_block = sysd.get("config") or {}
    ecology_cfg = None
    metabolism_cfg = None
    if cfg_block.get("ecology"):
        ecology_cfg = _resolve(config_dir, cfg_block["ecology"]) or Path(cfg_block["ecology"])
    if cfg_block.get("metabolism"):
        metabolism_cfg = _resolve(config_dir, cfg_block["metabolism"]) or Path(cfg_block["metabolism"])

    # Environments from registry
    env_specs = (sysd.get("registry") or {}).get("environments") or []
    env_files: List[Path] = []
    if env_specs:
        for item in env_specs:
            if isinstance(item, str):
                env_files.append(_resolve(envs_dir, item if item.endswith(".yml") else f"{item}.yml"))
            elif isinstance(item, dict):
                fid = item.get("file") or f"{item.get('id')}.yml"
                p = Path(fid)
                env_files.append(_resolve(envs_dir, fid) if not p.is_absolute() else p)
    else:
        env_files = sorted((envs_dir or root).glob("*.yml"))

    # Filter existing
    env_files = [p for p in env_files if p and p.exists()]

    return {
        "root": root,
        "system": sysd,
        "paths": paths,
        "ecology_cfg": ecology_cfg,
        "metabolism_cfg": metabolism_cfg,
        "environment_files": env_files,
    }

def iter_spot_files_for_env(env_file: Path, sys_paths: Dict[str, Any]) -> List[Tuple[str, Path]]:
    """
    Yield (spot_id, spot_path) for one environment file.
    Resolution rules:
      1) if environment.spots_dir is set, prefer that
      2) else if system.paths.spots_dir is set, use that
      3) else use a 'spots' subdir next to the env file
    If environment.spots is listed, resolve those; else glob *.yml in chosen dir.
    """
    env = _read_yaml(env_file).get("environment") or {}
    base = env_file.parent

    env_spots_dir = env.get("spots_dir")
    sys_spots_dir = sys_paths.get("spots_dir")
    if env_spots_dir:
        spots_base = _resolve(base, env_spots_dir)
    elif sys_spots_dir:
        spots_base = Path(sys_spots_dir)
    else:
        spots_base = base / "spots"

    out: List[Tuple[str, Path]] = []

    # Explicit list
    if env.get("spots"):
        for s in env["spots"]:
            if not isinstance(s, dict):
                continue
            sid = s.get("id") or s.get("name")
            f = s.get("file")
            if not sid or not f:
                continue
            p = Path(f)
            spath = (spots_base / p) if (p.parent == Path(".")) else _resolve(base, f)
            out.append((sid, spath))
        return out

    # Glob fallback
    if spots_base and spots_base.exists():
        for p in sorted(spots_base.glob("*.yml")):
            out.append((p.stem, p))
    return out

# ---------- Small helpers expected by CLI modules ----------

def read_spot_yaml(spot_path: Path) -> dict:
    """Load a spot YAML and return the dict under top-level 'spot' (or {})."""
    d = _read_yaml(Path(spot_path))
    return d.get("spot") or {}

def read_microbe_yaml(microbe_path: Path) -> dict:
    """Load a microbe YAML and return the dict under top-level 'microbe' (or {})."""
    d = _read_yaml(Path(microbe_path))
    return d.get("microbe") or {}

def load_microbe_registry(system_yml: Path) -> Dict[str, Path]:
    """
    From system.yml, return {microbe_id: microbe_yaml_path} with paths resolved
    via system.paths.microbes_dir (or system root).
    """
    sys_info = load_system(system_yml)
    root = sys_info["root"]
    sysd = sys_info["system"]
    paths = sys_info["paths"]
    microbes_dir = Path(paths.get("microbes_dir") or root / "microbes")

    reg = (sysd.get("registry") or {}).get("microbes") or []
    out: Dict[str, Path] = {}
    for item in reg:
        if isinstance(item, str):
            mid = item.replace(".yml", "")
            out[mid] = (microbes_dir / f"{mid}.yml")
        elif isinstance(item, dict) and item.get("id"):
            fid = item.get("file") or f"{item['id']}.yml"
            p = Path(fid)
            out[item["id"]] = (microbes_dir / p) if not p.is_absolute() else p
    # keep only existing files
    return {k: v for k, v in out.items() if v.exists()}
