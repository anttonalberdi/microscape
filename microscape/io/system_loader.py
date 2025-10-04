from __future__ import annotations
from pathlib import Path
from typing import List, Dict, Tuple
import yaml, glob

def _read_yaml(p: Path) -> dict:
    return yaml.safe_load(p.read_text())

def _resolve(base: Path, maybe_path: str | None) -> Path | None:
    if not maybe_path:
        return None
    p = Path(maybe_path)
    return (base / p) if not p.is_absolute() else p

def load_system(system_yml: Path) -> dict:
    system_yml = Path(system_yml).resolve()
    root = system_yml.parent
    sysd = _read_yaml(system_yml)["system"]

    paths = (sysd.get("paths") or {})
    env_dir = _resolve(root, paths.get("environments_dir")) or (root / "environments")
    # If there is an explicit registry list of env files, honor it; else glob.
    reg = sysd.get("registry") or {}
    env_list = reg.get("environments")
    env_files: List[Path] = []
    if env_list:
        for item in env_list:
            # allow either "E001" (infer environments/E001.yml) or a file path
            if isinstance(item, str) and item.endswith(".yml"):
                env_files.append(_resolve(root, item))
            else:
                env_files.append(env_dir / f"{item}.yml")
    else:
        env_files = sorted(Path(p) for p in glob.glob(str(env_dir / "*.yml")))

    # normalize and filter missing safely
    env_files = [p for p in env_files if p and p.exists()]

    return {
        "root": root,
        "system": sysd,
        "ecology_cfg": _resolve(root, (sysd.get("config") or {}).get("ecology")),
        "environment_files": env_files,
        "paths": paths,
    }

def iter_spot_files_for_env(env_file: Path, sys_paths: Dict) -> List[Tuple[str, Path]]:
    envd = _read_yaml(env_file)["environment"]
    base = env_file.parent
    spots_dir = envd.get("spots_dir")  # optional
    out: List[Tuple[str, Path]] = []
    if spots_dir:
        sdir = _resolve(base, spots_dir)
        for p in sorted(sdir.glob("*.yml")):
            sid = p.stem
            out.append((sid, p))
        return out
    # else read explicit list
    spots = envd.get("spots") or []
    for s in spots:
        sid = s.get("id")
        f = s.get("file")
        if sid and f:
            out.append((sid, _resolve(base, f)))
    return out
