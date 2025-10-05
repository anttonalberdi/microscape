# microscape/io/microbe_registry.py
from __future__ import annotations
from pathlib import Path
from typing import Dict, Tuple, List, Any
import yaml

def _read_yaml(p: Path) -> dict:
    with p.open("r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)

def build_microbe_model_map(system_root: Path, registry: dict) -> Tuple[Dict[str, Path], List[str]]:
    """
    Returns:
      microbe_models: {microbe_id -> absolute Path to SBML}
      warnings: list of warning strings
    """
    warnings: List[str] = []
    microbe_models: Dict[str, Path] = {}

    microbe_entries = (registry or {}).get("microbes") or []
    for entry in microbe_entries:
        if not isinstance(entry, dict) or "id" not in entry or "file" not in entry:
            warnings.append(f"Invalid microbe registry entry (need id,file): {entry}")
            continue

        mid = entry["id"]
        myml = (system_root / entry["file"]).resolve()
        if not myml.exists():
            warnings.append(f"Microbe YAML not found: {myml}")
            continue

        mdata = _read_yaml(myml)
        if not isinstance(mdata, dict) or "microbe" not in mdata:
            warnings.append(f"{myml}: missing top-level 'microbe' key.")
            continue

        m = mdata["microbe"]
        model = (m.get("model") or {})
        rel_model = model.get("path")
        if not rel_model:
            warnings.append(f"{myml}: model.path missing")
            continue

        sbml_path = (myml.parent / rel_model).resolve()
        if not sbml_path.exists():
            warnings.append(f"{myml}: model.path points to missing file: {rel_model}")
            continue

        microbe_models[mid] = sbml_path

    return microbe_models, warnings
