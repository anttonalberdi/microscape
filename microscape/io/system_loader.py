from __future__ import annotations
from pathlib import Path
from typing import List, Dict, Tuple, Any
import yaml

# ---------- helpers ----------

def _read_yaml(p: Path) -> dict:
    return yaml.safe_load(p.read_text())

def _resolve(base: Path, maybe: str | Path | None) -> Path | None:
    if not maybe:
        return None
    p = Path(maybe)
    return p if p.is_absolute() else (base / p)

# ---------- loaders ----------

def load_system(system_yml: Path) -> dict:
    """
    Load system.yml and resolve key paths relative to the directory containing system.yml.
    Returns:
      {
        "root": Path,                    # directory of system.yml
        "system": dict,                  # parsed 'system' section
        "paths": {                       # resolved absolute Paths
          "config_dir": Path,
          "environments_dir": Path,
          "spots_dir": Path,
          "microbes_dir": Path,
        },
        "ecology_cfg": Path|None,        # resolved ecology config (if provided)
        "environment_files": [Path],     # resolved environment YAMLs
      }
    """
    system_yml = Path(system_yml).resolve()
    root = system_yml.parent

    data = _read_yaml(system_yml)
    if not isinstance(data, dict) or "system" not in data:
        raise ValueError(f"{system_yml}: top-level 'system' key missing or invalid.")
    sysd: Dict[str, Any] = data["system"]

    # Resolve core directories (default to common names under root)
    paths_cfg: Dict[str, str] = sysd.get("paths") or {}
    config_dir     = _resolve(root, paths_cfg.get("config_dir"))       or (root / "config")
    envs_dir       = _resolve(root, paths_cfg.get("environments_dir")) or (root / "environments")
    spots_dir      = _resolve(root, paths_cfg.get("spots_dir"))        or (root / "spots")
    microbes_dir   = _resolve(root, paths_cfg.get("microbes_dir"))     or (root / "microbes")

    # Resolve ecology rules (relative to config_dir unless absolute)
    ecology_rel = (sysd.get("config") or {}).get("ecology")
    ecology_cfg = None
    if ecology_rel:
        ecology_cfg = _resolve(config_dir, ecology_rel)

    # Resolve environment files from registry (supports {id,file} or plain strings)
    env_specs = (sysd.get("registry") or {}).get("environments") or []
    env_files: List[Path] = []
    if env_specs:
        for item in env_specs:
            if isinstance(item, str):
                # "E001" or "E001.yml"
                f = f"{item}.yml" if not item.endswith(".yml") else item
                env_files.append((_resolve(envs_dir, f) or Path(f)).resolve())
            elif isinstance(item, dict):
                fid = item.get("file") or f"{item.get('id')}.yml"
                env_files.append((_resolve(envs_dir, fid) or Path(fid)).resolve())
    else:
        env_files = sorted(envs_dir.glob("*.yml"))

    # Keep only existing files
    env_files = [p for p in env_files if p.exists()]

    return {
        "root": root,
        "system": sysd,
        "paths": {
            "config_dir":   config_dir.resolve(),
            "environments_dir": envs_dir.resolve(),
            "spots_dir":    spots_dir.resolve(),
            "microbes_dir": microbes_dir.resolve(),
        },
        "ecology_cfg": ecology_cfg.resolve() if ecology_cfg else None,
        "environment_files": env_files,
    }

def iter_spot_files_for_env(env_file: Path, sys_paths: Dict[str, Path]) -> List[Tuple[str, Path]]:
    """
    List all spot files for a given environment file.

    IMPORTANT: Spot file paths are resolved relative to the **system-level** spots_dir,
    not the environment file's folder. This keeps all spots in a single project-level
    location and matches the design where system.yml is the single source of path truth.

    Returns: list of (spot_id, spot_path)
    """
    env_file = Path(env_file).resolve()
    env = _read_yaml(env_file).get("environment") or {}
    spots_base = Path(sys_paths.get("spots_dir") or (env_file.parent / "spots")).resolve()

    out: List[Tuple[str, Path]] = []

    # Explicit spot list in the environment
    spots = env.get("spots")
    if isinstance(spots, list) and spots:
        for s in spots:
            if not isinstance(s, dict):
                continue
            sid = s.get("id")
            f = s.get("file")
            if not sid or not f:
                continue
            p = Path(f)
            spath = p if p.is_absolute() else (spots_base / p)
            out.append((sid, spath.resolve()))
        return out

    # Otherwise, glob all *.yml under the system-level spots_dir
    if spots_base.exists():
        for p in sorted(spots_base.glob("*.yml")):
            out.append((p.stem, p.resolve()))
    return out
