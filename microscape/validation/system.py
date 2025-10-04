from __future__ import annotations
from pathlib import Path
from typing import Any, Dict, List, Tuple
import xml.etree.ElementTree as ET

try:
    import yaml  # PyYAML
except Exception:
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

def _parse_sbml_stats(model_file: Path, warnings: List[str]):
    """Return (species_total, unique_metabolite_ids_set, genes_set, expressed_genes_set)."""
    species_total = 0
    mets = set()
    genes = set()
    expressed = set()
    try:
        tree = ET.parse(model_file)
        root = tree.getroot()
        ns = {
            "sbml": "http://www.sbml.org/sbml/level3/version1/core",
            "fbc": "http://www.sbml.org/sbml/level3/version1/fbc/version2",
        }
        for sp in root.findall(".//sbml:listOfSpecies/sbml:species", ns):
            sid = sp.get("id") or ""
            species_total += 1
            base = sid.split("_")[0] if "_" in sid else sid
            if base:
                mets.add(base)
        for gp in root.findall(".//fbc:listOfGeneProducts/fbc:geneProduct", ns):
            gid = gp.get("{http://www.sbml.org/sbml/level3/version1/fbc/version2}id") or gp.get("id")
            if gid:
                genes.add(gid)
        for gpr in root.findall(".//fbc:geneProductAssociation", ns):
            for ref in gpr.findall(".//fbc:geneProductRef", ns):
                gid = ref.get("{http://www.sbml.org/sbml/level3/version1/fbc/version2}geneProduct")
                if gid:
                    expressed.add(gid)
    except Exception as e:
        warnings.append(f"{model_file}: SBML parse warning: {e}")
    return species_total, mets, genes, expressed

def validate_system(system_path: Path) -> Tuple[Dict[str, Any], List[str], List[str]]:
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

    env_entries = registry.get("environments", [])
    microbe_entries = registry.get("microbes", [])
    if not isinstance(env_entries, list):
        errors.append(f"{system_path}: system.registry.environments must be a list of {id,file} entries.")
        env_entries = []
    if not isinstance(microbe_entries, list):
        errors.append(f"{system_path}: system.registry.microbes must be a list of {id,file} entries.")
        microbe_entries = []

    def _dupes(ids: List[str]) -> List[str]:
        seen, d = set(), []
        for x in ids:
            if x in seen: d.append(x)
            seen.add(x)
        return d

    microbe_ids = [e.get("id") for e in microbe_entries if isinstance(e, dict)]
    env_ids = [e.get("id") for e in env_entries if isinstance(e, dict)]
    for name, ids in (("microbe", microbe_ids), ("environment", env_ids)):
        d = _dupes([x for x in ids if x is not None])
        if d:
            errors.append(f"{system_path}: duplicate {name} ids: {sorted(d)}")

    # Microbes
    microbe_index: Dict[str, Dict[str, str]] = {}
    model_files: List[Path] = []
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
            else:
                model_files.append(model_file)
        if not model_format:
            warnings.append(f"{myml}: model.format missing (expected 'sbml-l3v1-fbc2')")
        elif str(model_format).lower() != "sbml-l3v1-fbc2":
            warnings.append(f"{myml}: model.format '{model_format}' != 'sbml-l3v1-fbc2'")
        microbe_index[mid] = {"path": str(myml), "genus": tax.get("genus")}

    # Environments & Spots
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

    # Spots + spatial stats
    seen_spot_ids = set()
    has_positions = False
    spots_with_pos = 0
    z_values: List[float] = []

    def _spot_pos(spot_obj):
        pos = spot_obj.get("pos_um") if isinstance(spot_obj, dict) else None
        if isinstance(pos, (list, tuple)) and len(pos) >= 2:
            x, y = pos[0], pos[1]
            z = pos[2] if len(pos) >= 3 else None
            return x, y, z
        x = spot_obj.get("x_um") if isinstance(spot_obj, dict) else None
        y = spot_obj.get("y_um") if isinstance(spot_obj, dict) else None
        z = spot_obj.get("z_um") if isinstance(spot_obj, dict) else None
        if x is not None and y is not None:
            return x, y, z
        return None

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

        pos = _spot_pos(spot)
        if pos is not None:
            has_positions = True
            spots_with_pos += 1
            _, _, z = pos
            if z is not None:
                try:
                    z_values.append(float(z))
                except Exception:
                    pass

    # Aggregate SBML stats across microbes
    total_species = 0
    unique_mets = set()
    unique_genes = set()
    expressed_genes = set()
    for mf in model_files:
        s_cnt, mets, genes, expr = _parse_sbml_stats(mf, warnings)
        total_species += s_cnt
        unique_mets |= mets
        unique_genes |= genes
        expressed_genes |= expr

    # Dimensionality
    if not has_positions:
        dims = 0
        z_uniform = None
    else:
        if not z_values:
            dims = 2
            z_uniform = True
        else:
            zmin = min(z_values); zmax = max(z_values)
            z_uniform = (abs(zmax - zmin) < 1e-9)
            dims = 2 if z_uniform else 3

    summary = {
        "root": str(root),
        "microbes": {
            "count": len(microbe_entries),
            "genera": sorted({v.get("genus") for v in microbe_index.values() if v.get("genus")}),
        },
        "environments": {"count": len(env_entries)},
        "spots": {"count": len(seen_spot_ids)},
        "models": {
            "species_total": total_species,
            "metabolites_unique": len(unique_mets),
            "genes_total": len(unique_genes),
            "genes_expressed": len(expressed_genes),
            "genes_unused": max(0, len(unique_genes) - len(expressed_genes)),
        },
        "space": {
            "has_positions": has_positions,
            "spots_with_position": spots_with_pos,
            "dimensions": dims,
            "z_uniform": z_uniform,
        },
    }
    return summary, errors, warnings
