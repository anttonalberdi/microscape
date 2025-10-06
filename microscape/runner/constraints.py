# microscape/runner/constraints.py
from __future__ import annotations
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import csv
from dataclasses import dataclass

from cobra.io import read_sbml_model, write_sbml_model
from cobra import Model
from cobra.util.solver import linear_reaction_coefficients

from ..io.system_loader import load_system, iter_spot_files_for_env, read_spot_yaml, read_microbe_yaml
from ..io.metabolism_rules import (
    MetabolismRules,
    mM_to_uptake_bounds,
    reaction_scale_from_expr,
)

@dataclass
class ConstraintResult:
    spot_id: str
    microbe: str
    mode: str        # environmental | transcriptional | combined
    changed_ex: int
    changed_internal: int
    warnings: List[str]

def _apply_environmental(model: Model, spot_mets: Dict[str, float], rules: MetabolismRules) -> Tuple[int, List[str]]:
    changed = 0
    warns: List[str] = []
    fmap = rules.metabolite_map or {}
    for met_id, conc in (spot_mets or {}).items():
        ex_id = fmap.get(met_id)
        if not ex_id:
            continue
        rxn = model.reactions.get_by_id(ex_id) if ex_id in model.reactions else None
        if not rxn:
            warns.append(f"Missing EX reaction in model: {ex_id}")
            continue
        lb, ub = mM_to_uptake_bounds(float(conc), rules.uptake)
        rxn.lower_bound = lb
        rxn.upper_bound = ub
        changed += 1
    return changed, warns

def _apply_transcriptional(model: Model, gene_expr: Dict[str, float], rules: MetabolismRules) -> Tuple[int, List[str]]:
    changed = 0
    warns: List[str] = []
    for rxn in model.reactions:
        # Skip exchange reactions (id starts with EX_ in your toy schema)
        if rxn.id.startswith("EX_"):
            continue
        scale = reaction_scale_from_expr(rxn, gene_expr, rules.transcription)
        if scale >= 0.999:
            continue
        # scale both bounds symmetrically
        lb0, ub0 = rxn.lower_bound, rxn.upper_bound
        rxn.lower_bound = lb0 * scale
        rxn.upper_bound = ub0 * scale
        changed += 1
    return changed, warns

def _collect_spot_gene_expr(spot_obj: dict, microbe_id: str) -> Dict[str, float]:
    """
    Expecting in spot YAML:
      spot.measurements.transcripts.values.Mxxxx.Gyyyy: TPM
    """
    meas = (spot_obj.get("measurements") or {})
    tx = (meas.get("transcripts") or {})
    vals = (tx.get("values") or {})
    per_microbe = vals.get(microbe_id) or {}
    # keys are gene IDs; values numeric TPM
    return {str(k): float(v) for k, v in per_microbe.items()}

def constrain_one(
    system_yml: Path,
    rules: MetabolismRules,
    out_csv: Path,
    mode: str = "environmental",
    write_models: bool = False,
    models_outdir: Optional[Path] = None,
    verbose: bool = False
) -> List[ConstraintResult]:

    sys_info = load_system(system_yml)
    env_files = sys_info["environment_files"]

    # Build microbe → model path from registry
    microbe_models: Dict[str, Path] = {}
    for m in (sys_info["system"].get("registry") or {}).get("microbes", []):
        mid = m.get("id")
        mf = m.get("file")
        if not mid or not mf:
            continue
        microbes_dir = Path(sys_info["paths"]["microbes_dir"])
        myml = (microbes_dir / mf).resolve()
        md = read_microbe_yaml(myml)
        model_path = (myml.parent / md["microbe"]["model"]["path"]).resolve()
        md = read_microbe_yaml(myml)  # raises if bad
        model_path = (myml.parent / (md["microbe"]["model"]["path"])).resolve()
        microbe_models[mid] = model_path

    results: List[ConstraintResult] = []

    # CSV
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["spot_id", "microbe", "mode", "changed_ex", "changed_internal", "warnings"])

        for env_file in env_files:
            for sid, spot_path in iter_spot_files_for_env(env_file, sys_info["paths"]):
                spot = read_spot_yaml(spot_path)
                # spot metabolites
                mets = ((spot.get("measurements") or {}).get("metabolites") or {}).get("values") or {}
                for microbe, model_path in microbe_models.items():
                    try:
                        model = read_sbml_model(str(model_path))
                    except Exception as e:
                        w.writerow([sid, microbe, mode, 0, 0, f"SBML load error: {e}"])
                        continue

                    changed_ex = changed_internal = 0
                    warns_all: List[str] = []

                    if mode in ("environmental", "combined"):
                        ex_changed, warns = _apply_environmental(model, mets, rules)
                        changed_ex += ex_changed
                        warns_all += warns

                    if mode in ("transcriptional", "combined"):
                        gene_expr = _collect_spot_gene_expr(spot, microbe)
                        int_changed, warns = _apply_transcriptional(model, gene_expr, rules)
                        changed_internal += int_changed
                        warns_all += warns

                    # Optionally write constrained model
                    if write_models:
                        models_outdir = models_outdir or (out_csv.parent / "models_constrained")
                        models_outdir.mkdir(parents=True, exist_ok=True)
                        out_model = models_outdir / f"{sid}__{microbe}.xml"
                        try:
                            write_sbml_model(model, str(out_model))
                        except Exception as e:
                            warns_all.append(f"Write SBML failed: {e}")

                    # write row
                    warn_str = "; ".join(warns_all) if warns_all else ""
                    w.writerow([sid, microbe, mode, changed_ex, changed_internal, warn_str])

                    results.append(ConstraintResult(
                        spot_id=sid, microbe=microbe, mode=mode,
                        changed_ex=changed_ex, changed_internal=changed_internal,
                        warnings=warns_all
                    ))
    return results
