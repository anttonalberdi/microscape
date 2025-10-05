# microscape/io/microbe_registry.py
from __future__ import annotations
from pathlib import Path
from typing import Dict, List, Tuple, Any
import yaml

def _read_yaml(p: Path) -> dict:
    return yaml.safe_load(p.read_text())

def _resolve(base: Path, rel_or_abs: str) -> Path:
    p = Path(rel_or_abs)
    return p if p.is_absolute() else (base / p)

def build_microbe_model_map(
    root: Path,
    system_obj: dict,
    sys_paths: dict,
) -> Tuple[Dict[str, Path], List[str]]:
    """
    Returns:
      (model_map, warnings)
      model_map: { microbe_id -> SBML Path }
    """
    warnings: List[str] = []

    # Bases
    microbes_dir = sys_paths.get("microbes_dir")
    microbes_base = _resolve(root, microbes_dir) if microbes_dir else (root / "microbes")

    # Pull registry (if any)
    registry = (system_obj or {}).get("registry", {}) or {}
    microbe_entries = registry.get("microbes", []) or []

    # If registry empty, fallback to glob all YAMLs in microbes_base
    if not microbe_entries:
        microbe_files = sorted((microbes_base).glob("*.yml"))
        microbe_entries = [{"id": p.stem, "file": str(p.relative_to(root))} for p in microbe_files]

    model_map: Dict[str, Path] = {}

    for ent in microbe_entries:
        # Support "M0001" or "M0001.yml" (string) or {id:, file:}
        if isinstance(ent, str):
            if ent.endswith(".yml"):
                myml = _resolve(microbes_base, ent) if "/" not in ent else _resolve(root, ent)
                mid = Path(ent).stem
            else:
                myml = microbes_base / f"{ent}.yml"
                mid = ent
        elif isinstance(ent, dict):
            mid = ent.get("id")
            f = ent.get("file", f"{mid}.yml")
            # if f has a subdir, resolve against root; else against microbes_base
            myml = _resolve(root, f) if ("/" in f or "\\" in f) else (microbes_base / f)
        else:
            warnings.append(f"Invalid microbe registry entry: {ent}")
            continue

        if not myml.exists():
            warnings.append(f"Microbe YAML not found: {myml}")
            continue

        try:
            mdata = _read_yaml(myml)
        except Exception as e:
            warnings.append(f"YAML load failed for {myml}: {e}")
            continue

        m = (mdata or {}).get("microbe") or {}
        if not m:
            warnings.append(f"{myml}: missing top-level 'microbe' key.")
            continue

        # Read model path relative to the microbe YAML file
        model = (m.get("model") or {})
        model_path = model.get("path")
        if not model_path:
            warnings.append(f"{myml}: microbe.model.path missing")
            continue

        sbml = _resolve(myml.parent, model_path)
        if not sbml.exists():
            warnings.append(f"{myml}: model.path not found: {sbml}")
            continue

        # Use file name id if registry omitted id
        mid = mid or m.get("id") or myml.stem
        model_map[str(mid)] = sbml

    return model_map, warnings
