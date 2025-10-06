# microscape/io/system_loader.py
from __future__ import annotations
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import yaml

# ---------- utils ----------

def _read_yaml(p: Path) -> dict:
    with p.open("r", encoding="utf-8") as fh:
        return yaml.safe_load(fh) or {}

def _resolve(base: Path, maybe: Optional[str]) -> Optional[Path]:
    if not maybe:
        return None
    p = Path(maybe)
    return (base / p) if not p.is_absolute() else p

def _resolve_under_dir(root: Path, dir_from_paths: Optional[str], entry_file: str | None, fallback_dirname: str) -> Path:
    """
    Resolve a registry 'file' (e.g., 'M0001.yml') under the directory declared
    in system.paths (e.g. microbes_dir). Falls back to <root>/<fallback_dirname>.
    """
    base_dir = _resolve(root, dir_from_paths) or (root / fallback_dirname)
    if not entry_file:
        raise ValueError(f"Missing 'file' in registry entry for {fallback_dirname[:-1]}")
    p = Path(entry_file)
    return (base_dir / p) if not p.is_absolute() else p

# ---------- public API ----------

def load_system(system_yml: Path) -> dict:
    """
    Loads system.yml and returns:
      {
        "root": <Path>,
        "system": <dict>,
        "paths": <dict>,
        "ecology_cfg": <Path | None>,
        "metabolism_cfg": <Path | None>,
        "environment_files": [Path, ...],
        "microbe_files": [Path, ...],
      }
    All paths are resolved relative to system.yml and the configured dirs.
    """
    system_yml = Path(system_yml).resolve()
    root = system_yml.parent
    sysd = _read_yaml(system_yml).get("system") or {}

    paths: Dict = sysd.get("paths") or {}
    config_dir = _resolve(root, paths.get("config_dir")) or (root / "config")
    envs_dir   = _resolve(root, paths.get("environments_dir")) or (root / "environments")
    spots_dir  = _resolve(root, paths.get("spots_dir")) or (root / "spots")
    microbes_dir = _resolve(root, paths.get("microbes_dir")) or (root / "microbes")

    # Configs
    cfg_block = sysd.get("config") or {}
    ecology_cfg    = _resolve(config_dir, cfg_block.get("ecology")) if cfg_block.get("ecology") else None
    metabolism_cfg = _resolve(config_dir, cfg_block.get("metabolism")) if cfg_block.get("metabolism") else None

    # Registry
    reg = sysd.get("registry") or {}
    env_specs = reg.get("environments") or []
    mic_specs = reg.get("microbes") or []

    # Resolve environment files
    env_files: List[Path] = []
    for e in env_specs:
        if isinstance(e, str):
            # "E001" or "E001.yml"
            fname = e if e.endswith(".yml") else f"{e}.yml"
            env_files.append(envs_dir / fname)
        elif isinstance(e, dict):
            fid = e.get("file") or (f"{e.get('id')}.yml" if e.get("id") else None)
            if fid:
                env_files.append(_resolve_under_dir(root, paths.get("environments_dir"), fid, "environments"))
    env_files = [p.resolve() for p in env_files if p and p.exists()]

    # Resolve microbe files
    microbe_files: List[Path] = []
    for m in mic_specs:
        if isinstance(m, str):
            fname = m if m.endswith(".yml") else f"{m}.yml"
            microbe_files.append(microbes_dir / fname)
        elif isinstance(m, dict):
            fid = m.get("file") or (f"{m.get('id')}.yml" if m.get("id") else None)
            if fid:
                microbe_files.append(_resolve_under_dir(root, paths.get("microbes_dir"), fid, "microbes"))
    microbe_files = [p.resolve() for p in microbe_files if p and p.exists()]

    return {
        "root": root,
        "system": sysd,
        "paths": {
            "config_dir": str(config_dir),
            "environments_dir": str(envs_dir),
            "spots_dir": str(spots_dir),
            "microbes_dir": str(microbes_dir),
        },
        "ecology_cfg": ecology_cfg,
        "metabolism_cfg": metabolism_cfg,
        "environment_files": env_files,
        "microbe_files": microbe_files,
    }

def iter_spot_files_for_env(env_file: Path, sys_paths: Dict) -> List[Tuple[str, Path]]:
    """
    Returns [(spot_id, spot_path), ...] for an environment YAML.
    Resolution priority:
      1) env.spots_dir (relative to env file)
      2) system.paths.spots_dir (relative to system root)
      3) <env_file_dir>/spots
    If env lists explicit {id,file}, resolve those relative to chosen spots_dir.
    """
    env = _read_yaml(env_file).get("environment") or {}
    base = env_file.parent

    # Candidates
    env_spots_dir = env.get("spots_dir")
    sys_spots_dir = sys_paths.get("spots_dir")
    if env_spots_dir:
        spots_base = _resolve(base, env_spots_dir)
    elif sys_spots_dir:
        spots_base = Path(sys_spots_dir)
    else:
        spots_base = base / "spots"

    out: List[Tuple[str, Path]] = []
    if env.get("spots"):
        for s in env["spots"]:
            sid = s.get("id")
            f   = s.get("file")
            if not sid or not f:
                continue
            p = Path(f)
            spath = (spots_base / p) if not p.is_absolute() else p
            out.append((sid, spath.resolve()))
        return out

    if spots_base and spots_base.exists():
        for p in sorted(spots_base.glob("*.yml")):
            out.append((p.stem, p.resolve()))
    return out

# Convenience readers (used by CLI/runners)
def read_spot_yaml(spot_path: Path) -> dict:
    return _read_yaml(spot_path).get("spot") or {}

def read_environment_yaml(env_path: Path) -> dict:
    return _read_yaml(env_path).get("environment") or {}

def read_microbe_yaml(microbe_path: Path) -> dict:
    return _read_yaml(microbe_path)
