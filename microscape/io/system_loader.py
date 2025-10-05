# microscape/io/system_loader.py
from __future__ import annotations
from pathlib import Path
from typing import List, Dict, Tuple, Any
import yaml

def _read_yaml(p: Path) -> dict:
    return yaml.safe_load(p.read_text())

def _resolve(base: Path, maybe: str | None) -> Path | None:
    if not maybe:
        return None
    q = Path(maybe)
    return (base / q) if not q.is_absolute() else q

def load_system(system_yml: Path) -> dict:
    """
    Returns:
      {
        "root": Path,
        "system": dict,
        "paths": dict,
        "ecology_cfg": Optional[Path],
        "metabolism_cfg": Optional[Path],
        "environment_files": List[Path],
      }
    """
    system_yml = Path(system_yml).resolve()
    root = system_yml.parent
    sysd = _read_yaml(system_yml)["system"]

    paths: Dict[str, Any] = sysd.get("paths") or {}
    envs_dir = _resolve(root, paths.get("environments_dir")) or (root / "environments")
    config_dir = _resolve(root, paths.get("config_dir")) or (root / "config")

    # Resolve configs (relative to config_dir unless absolute)
    ecology_rel = (sysd.get("config") or {}).get("ecology")
    metabolism_rel = (sysd.get("config") or {}).get("metabolism")

    ecology_cfg = None
    if ecology_rel:
        ecology_cfg = _resolve(config_dir, ecology_rel) if not Path(ecology_rel).is_absolute() else Path(ecology_rel)

    metabolism_cfg = None
    if metabolism_rel:
        metabolism_cfg = _resolve(config_dir, metabolism_rel) if not Path(metabolism_rel).is_absolute() else Path(metabolism_rel)

    # Environments list
    env_specs = ((sysd.get("registry") or {}).get("environments") or [])
    env_files: List[Path] = []
    if env_specs:
        for item in env_specs:
            if isinstance(item, str):
                env_files.append(_resolve(envs_dir, item if item.endswith(".yml") else f"{item}.yml"))
            elif isinstance(item, dict):
                fid = item.get("file") or f"{item.get('id')}.yml"
                env_files.append(_resolve(envs_dir, fid) if not Path(fid).is_absolute() else Path(fid))
    else:
        env_files = sorted(envs_dir.glob("*.yml"))

    # Keep only existing files
    env_files = [p for p in env_files if p and p.exists()]

    return {
        "root": root,
        "system": sysd,
        "paths": paths,
        "ecology_cfg": ecology_cfg,
        "metabolism_cfg": metabolism_cfg,
        "environment_files": env_files,
    }

def iter_spot_files_for_env(env_file: Path, sys_paths: Dict) -> List[Tuple[str, Path]]:
    """
    Resolve spot YAML files for a given environment.

    Rules:
      - Prefer environment.spots list (id+file), relative to:
          * environment-level spots_dir if present
          * else system.paths.spots_dir if present
          * else default "spots" next to the environment YAML
      - If no explicit list, glob *.yml under the chosen spots_dir.
    """
    env = _read_yaml(env_file)["environment"]
    base = env_file.parent

    env_spots_dir = env.get("spots_dir")
    sys_spots_dir = sys_paths.get("spots_dir")

    if env_spots_dir:
        spots_base = _resolve(base, env_spots_dir)
    elif sys_spots_dir:
        # << matches your system.yml (paths.spots_dir is relative to system root).
        # For environments/E001.yml, spots live under <system_root>/<sys_spots_dir>.
        # base is <system_root>/environments, so resolve against system root:
        system_root = base.parent
        spots_base = _resolve(system_root, sys_spots_dir)
    else:
        spots_base = base / "spots"

    out: List[Tuple[str, Path]] = []

    if env.get("spots"):
        for s in env["spots"]:
            sid = s.get("id")
            f = s.get("file")
            if not sid or not f:
                continue
            q = Path(f)
            # If f is just a filename, prefix with spots_base.
            # If f has subdirs, still make it relative to environment base.
            spath = (spots_base / q) if q.parent == Path(".") else (base / q)
            out.append((sid, spath.resolve()))
        return out

    if spots_base and spots_base.exists():
        for p in sorted(spots_base.glob("*.yml")):
            out.append((p.stem, p.resolve()))
    return out
