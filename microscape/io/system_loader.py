from __future__ import annotations
from pathlib import Path
from typing import List, Dict, Tuple
import yaml, glob

def _read_yaml(p: Path) -> dict:
    return yaml.safe_load(p.read_text())

def _resolve(base: Path, maybe: str | None) -> Path | None:
    if not maybe:
        return None
    p = Path(maybe)
    return (base / p) if not p.is_absolute() else p

def load_system(system_yml: Path) -> dict:
    system_yml = Path(system_yml).resolve()
    root = system_yml.parent
    sysd = _read_yaml(system_yml)["system"]

    paths: Dict = sysd.get("paths") or {}
    envs_dir = _resolve(root, paths.get("environments_dir")) or (root / "environments")
    config_dir = _resolve(root, paths.get("config_dir")) or (root / "config")

    # Resolve ecology rules (e.g. "ecology.yml" under config_dir)
    ecology_rel = (sysd.get("config") or {}).get("ecology")
    ecology_cfg = None
    if ecology_rel:
        ecology_cfg = _resolve(config_dir, ecology_rel) if not Path(ecology_rel).is_absolute() else Path(ecology_rel)

    # Environments from registry (dicts or strings)
    env_specs = ((sysd.get("registry") or {}).get("environments") or [])
    env_files: List[Path] = []
    if env_specs:
        for item in env_specs:
            if isinstance(item, str):
                # "E001" or "E001.yml"
                if item.endswith(".yml"):
                    env_files.append(_resolve(envs_dir, item))
                else:
                    env_files.append(envs_dir / f"{item}.yml")
            elif isinstance(item, dict):
                # {id: E001, file: E001.yml}
                fid = item.get("file") or f"{item.get('id')}.yml"
                p = Path(fid)
                env_files.append(_resolve(envs_dir, fid) if not p.is_absolute() else p)
    else:
        env_files = sorted(envs_dir.glob("*.yml"))

    # filter to existing files
    env_files = [p for p in env_files if p and p.exists()]

    return {
        "root": root,
        "system": sysd,
        "paths": paths,
        "ecology_cfg": ecology_cfg,
        "environment_files": env_files,
    }

def iter_spot_files_for_env(env_file: Path, sys_paths: Dict) -> List[Tuple[str, Path]]:
    env = _read_yaml(env_file)["environment"]
    base = env_file.parent

    # Prefer per-environment spots_dir; else system-level; else "spots" next to env
    env_spots_dir = env.get("spots_dir")
    sys_spots_dir = sys_paths.get("spots_dir")
    if env_spots_dir:
        spots_base = _resolve(base, env_spots_dir)
    elif sys_spots_dir:
        spots_base = _resolve(base, sys_spots_dir)
    else:
        spots_base = base / "spots"

    out: List[Tuple[str, Path]] = []

    # If explicit list provided
    if env.get("spots"):
        for s in env["spots"]:
            sid = s.get("id")
            f = s.get("file")
            if not sid or not f:
                continue
            p = Path(f)
            # If only a filename, prefix with spots_base; if includes subdirs, still make it relative to base
            spath = (spots_base / p) if (p.parent == Path(".")) else _resolve(base, f)
            out.append((sid, spath))
        return out

    # Else, if a directory is specified, glob it
    if spots_base and spots_base.exists():
        for p in sorted(spots_base.glob("*.yml")):
            out.append((p.stem, p))
    return out
