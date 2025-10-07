# microscape/validate/system.py
from __future__ import annotations
from pathlib import Path
from typing import Any, Dict, List, Tuple, Set
import xml.etree.ElementTree as ET
import yaml

from ..io.system_loader import load_system, iter_spot_files_for_env

ALLOWED_MEASUREMENT_TYPES = {"counts", "relative", "biomass_gDW"}

def _read_yaml(p: Path) -> dict:
    with p.open("r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)

def _parse_sbml_stats(model_file: Path, warnings: List[str]):
    """Return (species_total, unique_metabolite_ids_set, genes_set, expressed_genes_set)."""
    species_total = 0
    mets: Set[str] = set()
    genes: Set[str] = set()
    expressed: Set[str] = set()
    try:
        tree = ET.parse(model_file)
        root = tree.getroot()
        ns = {
            "sbml": "http://www.sbml.org/sbml/level3/version1/core",
            "fbc": "http://www.sbml.org/sbml/level3/version1/fbc/version2",
        }
        # species
        for sp in root.findall(".//sbml:listOfSpecies/sbml:species", ns):
            sid = sp.get("id") or ""
            species_total += 1
            base = sid.split("_")[0] if "_" in sid else sid
            if base:
                mets.add(base)
        # genes
        for gp in root.findall(".//fbc:listOfGeneProducts/fbc:geneProduct", ns):
            gid = gp.get("{http://www.sbml.org/sbml/level3/version1/fbc/version2}id") or gp.get("id")
            if gid:
                genes.add(gid)
        # expressed in any GPR
        for gpr in root.findall(".//fbc:geneProductAssociation", ns):
            for ref in gpr.findall(".//fbc:geneProductRef", ns):
                gid = ref.get("{http://www.sbml.org/sbml/level3/version1/fbc/version2}geneProduct")
                if gid:
                    expressed.add(gid)
    except Exception as e:
        warnings.append(f"{model_file}: SBML parse warning: {e}")
    return species_total, mets, genes, expressed

def validate_system(system_yml: Path) -> Tuple[Dict[str, Any], List[str], List[str]]:
    """Validate a microscape project starting at system.yml."""
    errors: List[str] = []
    warnings: List[str] = []

    try:
        loaded = load_system(Path(system_yml))
    except Exception as e:
        return {}, [f"{system_yml}: {e}"], []

    root: Path               = loaded["root"]
    sysd: Dict[str, Any]     = loaded["system"]
    paths: Dict[str, Path]   = loaded["paths"]
    env_files: List[Path]    = loaded["environment_files"]
    mic_files: List[Path]    = loaded["microbe_files"]

    # Basic required keys in system section
    for k in ["id", "type", "registry"]:
        if k not in sysd:
            errors.append(f"{system_yml}: missing system.{k}")

    # ---- validate microbe registry files ----
    microbe_index: Dict[str, Dict[str, Any]] = {}
    model_files: List[Path] = []
    for f in mic_files:
        try:
            mdata = _read_yaml(f)
        except FileNotFoundError:
            errors.append(f"File not found: {f}")
            continue
        except Exception as e:
            errors.append(f"YAML parse error in {f}: {e}")
            continue

        if not isinstance(mdata, dict) or "microbe" not in mdata:
            errors.append(f"{f}: missing top-level 'microbe' key.")
            continue
        m = mdata["microbe"]
        mid = m.get("id") or f.stem
        tax = m.get("taxonomy", {}) or {}
        model = (m.get("model") or {})

        # taxonomy warnings (non-fatal)
        for k in ["phylum","class","order","family","genus","species"]:
            if k not in tax:
                warnings.append(f"{f}: taxonomy.{k} missing")

        # SBML model presence
        model_path = model.get("path")
        if not model_path:
            warnings.append(f"{f}: model.path missing (OK if not using SBML yet)")
        else:
            mf = (f.parent / model_path) if not Path(model_path).is_absolute() else Path(model_path)
            if mf.exists():
                model_files.append(mf)
            else:
                errors.append(f"{f}: model.path points to missing file: {model_path}")

        microbe_index[mid] = {"path": str(f), "genus": tax.get("genus")}

    # ---- validate environments & spots (resolve via system-level spots_dir) ----
    spot_files: List[Path] = []
    seen_spot_ids: Set[str] = set()
    has_positions = False
    spots_with_pos = 0
    z_values: List[float] = []

    def _spot_pos(spot_obj: Dict[str, Any]):

        def extract(obj: Dict[str, Any] | None):
            if not isinstance(obj, dict):
                return None
            # list-style
            pos = obj.get("pos_um") or obj.get("position_um") or obj.get("pos")
            if isinstance(pos, (list, tuple)) and len(pos) >= 2:
                x, y = pos[0], pos[1]
                z = pos[2] if len(pos) >= 3 else None
                return x, y, z
            # key-style
            x = obj.get("x_um", obj.get("x"))
            y = obj.get("y_um", obj.get("y"))
            z = obj.get("z_um", obj.get("z"))
            if x is not None and y is not None:
                return x, y, z
            return None

        # Prefer nested 'position' if present, then fall back to top-level
        nested = extract(spot_obj.get("position"))
        if nested is not None:
            return nested
        return extract(spot_obj)

    for ef in env_files:
        try:
            edata = _read_yaml(ef)
        except FileNotFoundError:
            errors.append(f"File not found: {ef}")
            continue
        except Exception as e:
            errors.append(f"YAML parse error in {ef}: {e}")
            continue

        if not isinstance(edata, dict) or "environment" not in edata:
            errors.append(f"{ef}: missing top-level 'environment' key.")
            continue
        env = edata["environment"]

        # Collect spot files for this environment
        for sid, sf in iter_spot_files_for_env(ef, paths):
            if not sf.exists():
                errors.append(f"{ef}: spot file not found: {sf}")
                continue
            spot_files.append(sf)

            # Read + inspect spot file
            try:
                sdata = _read_yaml(sf)
            except Exception as e:
                errors.append(f"YAML parse error in {sf}: {e}")
                continue
            if not isinstance(sdata, dict) or "spot" not in sdata:
                errors.append(f"{sf}: missing top-level 'spot' key.")
                continue
            spot = sdata["spot"]

            # ID uniqueness
            spot_id = spot.get("name") or spot.get("id") or sid
            if spot_id in seen_spot_ids:
                warnings.append(f"{sf}: duplicate spot id/name: {spot_id}")
            seen_spot_ids.add(spot_id)

            # Position
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

            # Measurements: microbes
            meas = spot.get("measurements", {})
            if not isinstance(meas, dict):
                warnings.append(f"{sf}: 'measurements' is not a mapping")
                continue

            microbes = meas.get("microbes")
            if not microbes:
                warnings.append(f"{sf}: no microbes measurements present")
            else:
                mtype = microbes.get("type")
                if mtype not in ALLOWED_MEASUREMENT_TYPES:
                    warnings.append(f"{sf}: microbes.type '{mtype}' not in {sorted(ALLOWED_MEASUREMENT_TYPES)}")
                vals = microbes.get("values", {}) or {}
                for mid in list(vals.keys()):
                    if mid not in microbe_index:
                        warnings.append(f"{sf}: microbes.values contains unknown microbe id '{mid}' (not in registry)")

    # ---- Aggregate SBML stats across microbes ----
    total_species = 0
    unique_mets: Set[str] = set()
    unique_genes: Set[str] = set()
    expressed_genes: Set[str] = set()
    for mf in model_files:
        s_cnt, mets, genes, expr = _parse_sbml_stats(mf, warnings)
        total_species += s_cnt
        unique_mets |= mets
        unique_genes |= genes
        expressed_genes |= expr

    # ---- Spatial dimensionality ----
    if not has_positions:
        dims = 0
        z_uniform = None
    else:
        if not z_values:
            dims = 2
            z_uniform = True
        else:
            zmin, zmax = min(z_values), max(z_values)
            z_uniform = (abs(zmax - zmin) < 1e-9)
            dims = 2 if z_uniform else 3

    summary = {
        "root": str(root),
        "microbes": {
            "count": len(microbe_index),
            "genera": sorted({v.get("genus") for v in microbe_index.values() if v.get("genus")}),
        },
        "environments": {"count": len(env_files)},
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
        "paths": {k: str(v) for k, v in paths.items()},
    }

    return summary, errors, warnings
