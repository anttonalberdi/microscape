from __future__ import annotations
from pathlib import Path
from typing import Any, Dict, List, Tuple

try:
    import yaml  # PyYAML
except Exception as e:
    # Delay import error until used by CLI to keep module importable in docs builds
    yaml = None  # type: ignore

ALLOWED_MEASUREMENT_TYPES = {"counts", "relative", "biomass_gDW"}

def _ensure_yaml():
    if yaml is None:
        raise RuntimeError("PyYAML is required. Install with: pip install pyyaml")

def _load_yaml(p: Path, errors: List[str]) -> Any:
    _ensure_yaml()
    try:
        with p.open("r", encoding="utf-8") as fh:
            return yaml.safe_load(fh)
    except FileNotFoundError:
        errors.append(f"File not found: {p}")
    except Exception as e:
        errors.append(f"YAML parse error in {p}: {e}")
    return None

def validate_system(system_path: Path) -> Tuple[Dict[str, Any], List[str], List[str]]:
    """
    Validate a microscape project starting from system.yml.

    Returns:
        summary: {root, microbes:{count,genera}, environments:{count}, spots:{count}}
        errors:  list of fatal issues (non-zero exit in CLI)
        warnings:list of non-fatal issues
    """
    errors: List[str] = []
    warnings: List[str] = []
    root = system_path.parent

    data = _load_yaml(system_path, errors)
    if not isinstance(data, dict) or "system" not in data:
        errors.append(f"{system_path}: top-level 'system' key missing or invalid.")
        return {}, errors, warnings
    sysobj = data["system"]

    # Required root keys
    for k in ["id", "type", "design", "registry"]:
        if k not in sysobj:
            errors.append(f"{system_path}: missing system.{k}")

    registry = sysobj.get("registry", {})
    if not isinstance(registry, dict):
        errors.append(f"{system_path}: system.registry must be a mapping.")
        registry = {}

    # Collect registry sections
    env_entries = registry.get("environments", [])
    microbe_entries = registry.get("microbes", [])
    if not isinstance(env_entries, list):
        errors.append(f"{system_path}: system.registry.environments must be a list of {{id,file}} entries.")
        env_entries = []
    if not isinstance(microbe_entries, list):
        errors.append(f"{system_path}: system.registry.microbes must be a list of {{id,file}} entries.")
        microbe_entries = []

    # Duplicate ID checks
    def _dupes(ids: List[str]) -> List[str]:
        seen, d = set(), []
        for x in ids:
            if x in seen:
                d.append(x)
            seen.add(x)
        return d

    microbe_ids = [e.get("id") for e in microbe_entries if isinstance(e, dict)]
    env_ids = [e.get("id") for e in env_entries if isinstance(e, dict)]
    for name, ids in (("microbe", microbe_ids), ("environment", env_ids)):
        d = _dupes([x for x in ids if x is not None])
        if d:
            errors.append(f"{system_path}: duplicate {name} ids: {sorted(d)}")

    # Validate Microbes
    microbe_index: Dict[str, Dict[str, str]] = {}
    for e in microbe_entries:
        if not isinstance(e, dict) or "id" not in e or "file" not in e:
            errors.append(f"{system_path}: invalid microbe registry entry (need id,file): {e}")
            continue
        mid = e["id"]
        myml = root / e["file"]
        mdata = _load_yaml(myml, errors)
        if not mdata or "microbe" not in mdata:
            errors.append(f"{myml}: missing top-level 'microbe' key.")
            continue
        m = mdata["microbe"]
        for k in ["id", "taxonomy", "model"]:
            if k not in m:
                errors.append(f"{myml}: missing microbe.{k}")
        if m.get("id") and m["id"] != mid:
            warnings.append(f"{myml}: id mismatch. registry={mid} file={m.get('id')}")
        tax = m.get("taxonomy", {}) or {}
        for k in ["phylum","class","order","family","genus","species"]:
            if k not in tax:
                warnings.append(f"{myml}: taxonomy.{k} missing")
        model = (m.get("model") or {})
        model_path = model.get("path")
        model_format = model.get("format")
        if not model_path:
            errors.append(f"{myml}: model.path missing")
        else:
            model_file = (myml.parent / model_path)
            if not model_file.exists():
                errors.append(f"{myml}: model.path points to missing file: {model_path}")
        if not model_format:
            warnings.append(f"{myml}: model.format missing (expected 'sbml-l3v1-fbc2')")
        elif str(model_format).lower() != "sbml-l3v1-fbc2":
            warnings.append(f"{myml}: model.format '{model_format}' != 'sbml-l3v1-fbc2'")
        microbe_index[mid] = {"path": str(myml), "genus": tax.get("genus")}

    # Validate Environments & Spots
    spot_files: List[Path] = []
    for e in env_entries:
        if not isinstance(e, dict) or "id" not in e or "file" not in e:
            errors.append(f"{system_path}: invalid environment registry entry (need id,file): {e}")
            continue
        eid = e["id"]
        eyml = root / e["file"]
        edata = _load_yaml(eyml, errors)
        if not edata or "environment" not in edata:
            errors.append(f"{eyml}: missing top-level 'environment' key.")
            continue
        env = edata["environment"]
        if env.get("id") and env["id"] != eid:
            warnings.append(f"{eyml}: id mismatch. registry={eid} file={env.get('id')}")
        spots = env.get("spots", [])
        if not isinstance(spots, list) or not spots:
            warnings.append(f"{eyml}: no spots listed under environment.spots")
            spots = []
        for sp in spots:
            if not isinstance(sp, dict) or "id" not in sp or "file" not in sp:
                errors.append(f"{eyml}: invalid spot entry (need id,file): {sp}")
                continue
            sf = root / sp["file"]
            if not sf.exists():
                errors.append(f"{eyml}: spot file not found: {sp['file']}")
            else:
                spot_files.append(sf)

    # Validate Spots
    seen_spot_ids = set()
    for sf in spot_files:
        sdata = _load_yaml(sf, errors)
        if not sdata or "spot" not in sdata:
            errors.append(f"{sf}: missing top-level 'spot' key.")
            continue
        spot = sdata["spot"]
        sid = spot.get("name") or spot.get("id")
        if not sid:
            warnings.append(f"{sf}: spot has no 'name' or 'id'")
        if sid in seen_spot_ids:
            warnings.append(f"{sf}: duplicate spot id/name: {sid}")
        seen_spot_ids.add(sid)
        meas = spot.get("measurements", {})
        if not isinstance(meas, dict):
            warnings.append(f"{sf}: measurements not a mapping")
        microbes = meas.get("microbes")
        if microbes:
            mtype = microbes.get("type")
            if mtype not in ALLOWED_MEASUREMENT_TYPES:
                warnings.append(f"{sf}: microbes.type '{mtype}' not in {sorted(ALLOWED_MEASUREMENT_TYPES)}")
            vals = microbes.get("values", {}) or {}
            for mid in list(vals.keys()):
                if mid not in microbe_index:
                    warnings.append(f"{sf}: microbes.values contains unknown microbe id '{mid}' (not in system.registry.microbes)")
        else:
            warnings.append(f"{sf}: no microbes measurements present")

    summary = {
        "root": str(root),
        "microbes": {
            "count": len(microbe_entries),
            "genera": sorted({v.get("genus") for v in microbe_index.values() if v.get("genus")}),
        },
        "environments": {"count": len(env_entries)},
        "spots": {"count": len(seen_spot_ids)},
    }
    return summary, errors, warnings
