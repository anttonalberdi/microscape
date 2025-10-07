# microscape/validate/system.py
from __future__ import annotations
from pathlib import Path
from typing import Any, Dict, List, Tuple
import xml.etree.ElementTree as ET
from ..io.system_loader import load_system, iter_spot_files_for_env, read_spot_yaml, read_microbe_yaml

def _parse_sbml_stats(model_file: Path, warnings: List[str]):
    """Return (species_total, unique_metabolite_bases_set, genes_set, expressed_genes_set)."""
    species_total = 0
    mets = set()
    genes = set()
    expressed = set()
    try:
        tree = ET.parse(model_file)
        root = tree.getroot()
        ns = {
            "sbml": "http://www.sbml.org/sbml/level3/version1/core",
            "fbc":  "http://www.sbml.org/sbml/level3/version1/fbc/version2",
        }
        for sp in root.findall(".//sbml:listOfSpecies/sbml:species", ns):
            sid = sp.get("id") or ""
            species_total += 1
            base = sid.split("_")[0] if "_" in sid else sid
            if base:
                mets.add(base)
        for gp in root.findall(".//fbc:listOfGeneProducts/fbc:geneProduct", ns):
            gid = gp.get("{http://www.sbml.org/sbml/level3/version1/fbc/version2}id") or gp.get("id")
            if gid: genes.add(gid)
        for gpr in root.findall(".//fbc:geneProductAssociation", ns):
            for ref in gpr.findall(".//fbc:geneProductRef", ns):
                gid = ref.get("{http://www.sbml.org/sbml/level3/version1/fbc/version2}geneProduct")
                if gid: expressed.add(gid)
    except Exception as e:
        warnings.append(f"{model_file}: SBML parse warning: {e}")
    return species_total, mets, genes, expressed

def validate_system(system_path: Path) -> Tuple[Dict[str, Any], List[str], List[str]]:
    errors: List[str] = []
    warnings: List[str] = []

    sys_info = load_system(Path(system_path))
    root = Path(sys_info["root"])
    env_files: List[Path] = sys_info["environment_files"]
    microbe_registry: Dict[str, Path] = sys_info["microbe_registry"]

    # Count microbes (registered), and accumulate SBML stats
    total_species = 0
    unique_mets = set()
    unique_genes = set()
    expressed_genes = set()

    microbe_count = len(microbe_registry)
    for mid, myml in microbe_registry.items():
        md = read_microbe_yaml(mid, sys_info)
        if not md or "microbe" not in md:
            errors.append(f"{myml}: missing top-level 'microbe' key.")
            continue
        model = (md["microbe"].get("model") or {})
        model_path = model.get("path")
        if not model_path:
            errors.append(f"{myml}: model.path missing")
            continue
        mp = (Path(myml).parent / model_path).resolve()
        if not mp.exists():
            errors.append(f"{myml}: model.path points to missing file: {model_path}")
            continue
        s_cnt, mets, genes, expr = _parse_sbml_stats(mp, warnings)
        total_species += s_cnt
        unique_mets |= mets
        unique_genes |= genes
        expressed_genes |= expr

    # Environments + spots
    seen_spot_ids = set()
    has_positions = False
    spots_with_pos = 0
    z_vals: List[float] = []

    for env_file in env_files:
        env = (read_spot_yaml(env_file) or {}).get("environment")
        if not env:
            errors.append(f"{env_file}: missing top-level 'environment' key.")
            continue

        for sid, spath in iter_spot_files_for_env(env_file, sys_info["paths"]):
            if not spath.exists():
                errors.append(f"{env_file}: spot file not found: {spath}")
                continue
            spot = read_spot_yaml(spath).get("spot", {})
            sid2 = spot.get("name") or spot.get("id") or sid
            if sid2 in seen_spot_ids:
                warnings.append(f"{spath}: duplicate spot id/name: {sid2}")
            seen_spot_ids.add(sid2)

            # microbes existence check
            microbes = ((spot.get("measurements") or {}).get("microbes") or {})
            vals = microbes.get("values", {}) or {}
            for mid in list(vals.keys()):
                if mid not in microbe_registry:
                    warnings.append(f"{spath}: microbes.values contains unknown microbe id '{mid}' (not in registry)")

            # positions
            pos = spot.get("position") or spot.get("pos_um")
            if isinstance(pos, dict):
                x = pos.get("x"); y = pos.get("y"); z = pos.get("z", None)
                if x is not None and y is not None:
                    has_positions = True; spots_with_pos += 1
                    if z is not None:
                        try: z_vals.append(float(z))
                        except: pass

    # dims
    if not has_positions:
        dims = 0; z_uniform = None
    else:
        if not z_vals:
            dims = 2; z_uniform = True
        else:
            zmin, zmax = min(z_vals), max(z_vals)
            z_uniform = abs(zmax - zmin) < 1e-9
            dims = 2 if z_uniform else 3

    summary = {
        "root": str(root),
        "microbes": {"count": microbe_count},
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
    }
    return summary, errors, warnings
