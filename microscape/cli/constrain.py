# microscape/cli/constrain.py
from __future__ import annotations
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import json, csv, typer
from rich.progress import Progress

from ..io.system_loader import load_system, iter_spot_files_for_env, read_spot_yaml, read_microbe_yaml
from ..io.metabolism_rules import load_rules, MetabolismRules

app = typer.Typer(add_completion=False, no_args_is_help=True)

# ----------------------- local helpers (avoid importing from loader) -----------------------

def _spot_metabolites(spot: dict) -> Dict[str, float]:
    """Return {metabolite_id: concentration} from a spot YAML dict."""
    meas = (spot or {}).get("measurements") or {}
    mets = (meas.get("metabolites") or {}).get("values") or {}
    # normalize values to float
    out: Dict[str, float] = {}
    for k, v in mets.items():
        try:
            out[str(k)] = float(v)
        except Exception:
            pass
    return out

def _spot_transcripts_by_microbe(spot: dict) -> Dict[str, Dict[str, float]]:
    """
    Return {microbe_id: {gene_id: expr_value}} from a spot YAML dict.
    Accepts either:
      transcripts:
        values:
          M0001: { G0001: 10, ... }
          M0002: ...
    """
    meas = (spot or {}).get("measurements") or {}
    tx_block = meas.get("transcripts") or {}
    values = tx_block.get("values") or {}
    out: Dict[str, Dict[str, float]] = {}
    for mid, genes in values.items():
        gmap: Dict[str, float] = {}
        if isinstance(genes, dict):
            for gid, val in genes.items():
                try:
                    gmap[str(gid)] = float(val)
                except Exception:
                    pass
        out[str(mid)] = gmap
    return out

def _rxn_kind(r) -> str:
    """Classify reaction kind: 'exchange' if EX_* or has any boundary metabolite; else 'internal'."""
    if r.id.startswith("EX_"):
        return "exchange"
    try:
        # cobra 0.29: metabolite has 'boundary' attribute in some loaders;
        # safer: check if any species id ends with '_e' and is boundaryCondition in SBML?
        # Here, approximate: if any metabolite is in extracellular compartment 'e', we treat EX_* explicitly above.
        pass
    except Exception:
        pass
    return "internal"

def _base_bounds(model) -> Dict[str, Tuple[float, float, str]]:
    out = {}
    for r in model.reactions:
        out[r.id] = (float(r.lower_bound), float(r.upper_bound), _rxn_kind(r))
    return out

# ----------------------- apply constraints -----------------------

def _apply_environmental_bounds(model, env_bounds: Dict[str, Tuple[float,float]]) -> List[str]:
    changed = []
    for rxn_id, (lb, ub) in env_bounds.items():
        r = model.reactions.get_by_id(rxn_id) if rxn_id in model.reactions else None
        if r is None:
            continue
        if (r.lower_bound != lb) or (r.upper_bound != ub):
            r.lower_bound = lb
            r.upper_bound = ub
            changed.append(rxn_id)
    return changed

def _apply_transcriptional_scaling(model, gexpr: Dict[str, float], mode: str, min_scale: float,
                                   norm: str, percentile: float) -> Tuple[List[str], List[str]]:
    """
    Scale INTERNAL reaction bounds from gene expression (never exchanges).
    Returns (changed_internal, warnings).
    """
    changed_internal: List[str] = []
    warnings: List[str] = []

    # Prepare gene score lookup (simple identity; model.genes are Gene objects with .id)
    # E-Flux: scale upper bounds by normalized gene score; lower bounds for reversible reactions scaled symmetrically.
    # Threshold: if expr < threshold, set ub=0 (or to min_scale*ub0).
    # We keep EX_* untouched.
    import numpy as np

    if not model.genes:
        return changed_internal, warnings

    # Collect expression array for genes present in the model
    vals = []
    for g in model.genes:
        x = gexpr.get(g.id, 0.0)
        try:
            vals.append(float(x))
        except Exception:
            vals.append(0.0)
    arr = np.array(vals, dtype=float)

    # Normalize
    if mode == "eflux":
        if norm == "max":
            denom = float(arr.max()) if arr.size and arr.max() > 0 else 1.0
        elif norm == "percentile":
            denom = float(np.percentile(arr, percentile)) if arr.size else 1.0
            if denom <= 0: denom = 1.0
        else:
            denom = float(arr.max()) if arr.size and arr.max() > 0 else 1.0
        # gene score: expr / denom, clipped to [min_scale, 1]
        gene_score = {g.id: max(min_scale, float(gexpr.get(g.id, 0.0)) / denom) for g in model.genes}
        def score_for_rxn(r):
            # if no GPR, skip
            if r.gene_reaction_rule.strip() == "":
                return None
            # simple OR-of-genes max score
            s = 0.0
            for g in r.genes:
                s = max(s, gene_score.get(g.id, min_scale))
            return max(min_scale, min(1.0, s))
        for r in model.reactions:
            if _rxn_kind(r) == "exchange":
                continue
            s = score_for_rxn(r)
            if s is None:
                continue
            lb0, ub0 = float(r.lower_bound), float(r.upper_bound)
            # scale only positive capacity; keep sign for reversible
            new_lb = lb0 * s if lb0 < 0 else (-abs(lb0) * s if lb0 < 0 else lb0 * s if lb0 > 0 else lb0*s)
            new_ub = ub0 * s if ub0 > 0 else ub0 * s
            # cleaner:
            if lb0 < 0: new_lb = lb0 * s
            if ub0 > 0: new_ub = ub0 * s
            # apply
            if (new_lb != r.lower_bound) or (new_ub != r.upper_bound):
                r.lower_bound = new_lb
                r.upper_bound = new_ub
                changed_internal.append(r.id)

    elif mode == "threshold":
        thr = float(norm)  # here we pass threshold through 'norm' param for simplicity
        for r in model.reactions:
            if _rxn_kind(r) == "exchange":
                continue
            if r.gene_reaction_rule.strip() == "":
                continue
            # if all genes for this rxn below threshold, clamp to min_scale*original
            all_below = True
            for g in r.genes:
                if float(gexpr.get(g.id, 0.0)) >= thr:
                    all_below = False
                    break
            if all_below:
                lb0, ub0 = float(r.lower_bound), float(r.upper_bound)
                new_lb = lb0 * min_scale if lb0 < 0 else 0.0
                new_ub = ub0 * min_scale if ub0 > 0 else 0.0
                if (new_lb != r.lower_bound) or (new_ub != r.upper_bound):
                    r.lower_bound = new_lb
                    r.upper_bound = new_ub
                    changed_internal.append(r.id)
    else:
        warnings.append(f"Unknown transcriptional mode '{mode}'")

    return changed_internal, warnings

# ----------------------- CLI command -----------------------

@app.command("constrain")
def constrain_cmd(
    system_yml: Path = typer.Argument(..., help="Path to system.yml"),
    outdir: Path = typer.Option("outputs/constrain", help="Output directory"),
    mode: str = typer.Option("environmental", help="environmental | transcriptional | combined"),
    verbose: bool = typer.Option(False, "--verbose", "-v"),
):
    """
    Build per-spot×microbe constraint proposals from:
      - environmental: metabolite concentrations → EX bounds (via metabolism rules)
      - transcriptional: gene expression → scale INTERNAL reaction capacity (never EX)
      - combined: apply both in one pass (EX from environment, INTERNAL from transcripts)
    Writes:
      - constraints__{mode}.tsv  (summary rows)
      - constraints__{mode}__reactions.json  (full per-reaction bounds for each spot×microbe)
    """
    sys_info = load_system(system_yml)
    metab_cfg = sys_info.get("metabolism_cfg")
    if not metab_cfg or not Path(metab_cfg).exists():
        typer.secho("Metabolism rules not found. Ensure system.config.metabolism points to a file.", fg=typer.colors.RED)
        raise typer.Exit(1)
    rules: MetabolismRules = load_rules(Path(metab_cfg))

    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    env_files: List[Path] = sys_info["environment_files"]

    summary_rows: List[Dict] = []
    reactions_payload: List[Dict] = []

    # progress bar
    total_spots = sum(len(iter_spot_files_for_env(e, sys_info["paths"])) for e in env_files) or 1
    typer.echo(f"Constraining (mode={mode})…")
    with Progress() as prog:
        task = prog.add_task("[cyan]Constraining…", total=total_spots)

        for env_file in env_files:
            for sid, spath in iter_spot_files_for_env(env_file, sys_info["paths"]):
                spot = read_spot_yaml(spath).get("spot", {})

                # gather inputs
                met_vals = _spot_metabolites(spot)            # {C0001: 5.2, ...}
                tx_vals  = _spot_transcripts_by_microbe(spot) # {M0001: {G0001: TPM, ...}, ...}

                # microbes in this spot (if absent, skip)
                microbes_vals = (((spot.get("measurements") or {}).get("microbes") or {}).get("values") or {})
                for mid in microbes_vals.keys():
                    myml = read_microbe_yaml(mid, sys_info)
                    if not myml or "microbe" not in myml:
                        if verbose: typer.echo(f"WARN: Microbe YAML not found for {mid}")
                        continue
                    model_rel = (myml["microbe"].get("model") or {}).get("path")
                    if not model_rel:
                        if verbose: typer.echo(f"WARN: No model.path in microbe {mid}")
                        continue
                    model_path = (Path(sys_info["microbe_registry"][mid]).parent / model_rel).resolve()
                    if not model_path.exists():
                        if verbose: typer.echo(f"WARN: model file missing for {mid}: {model_path}")
                        continue

                    import cobra
                    model = cobra.io.read_sbml_model(str(model_path))

                    # collect original bounds for save
                    base = _base_bounds(model)

                    changed_ex: List[str] = []
                    changed_int: List[str] = []
                    warns: List[str] = []

                    # apply per mode
                    if mode in ("environmental", "combined"):
                        env_bounds: Dict[str, Tuple[float,float]] = {}
                        for met_id, conc in met_vals.items():
                            ex_id = rules.metabolite_map.get(met_id)
                            if not ex_id:
                                continue
                            lb_env, ub_env = rules.uptake_to_bounds(conc)
                            env_bounds[ex_id] = (lb_env, ub_env)
                        changed_ex = _apply_environmental_bounds(model, env_bounds)

                    if mode in ("transcriptional", "combined"):
                        gexpr = tx_vals.get(mid, {})
                        tmode = (rules.transcription.mode or "eflux")
                        if tmode == "eflux":
                            cint, w = _apply_transcriptional_scaling(
                                model, gexpr,
                                mode="eflux",
                                min_scale=rules.transcription.min_scale,
                                norm=rules.transcription.eflux_norm,
                                percentile=rules.transcription.percentile,
                            )
                            changed_int.extend(cint); warns.extend(w)
                        elif tmode == "threshold":
                            cint, w = _apply_transcriptional_scaling(
                                model, gexpr,
                                mode="threshold",
                                min_scale=rules.transcription.min_scale,
                                norm=str(rules.transcription.threshold_TPM),
                                percentile=0.0,
                            )
                            changed_int.extend(cint); warns.extend(w)
                        else:
                            warns.append(f"Unknown transcription mode: {tmode}")

                    # summary row (one per spot×microbe)
                    summary_rows.append({
                        "spot_id": sid,
                        "microbe": mid,
                        "mode": mode,
                        "changed_ex": len(changed_ex),
                        "changed_internal": len(changed_int),
                        "warnings": "; ".join(warns) if warns else "",
                    })

                    # full reaction bounds payload
                    after = _base_bounds(model)
                    reactions_payload.append({
                        "spot_id": sid,
                        "microbe": mid,
                        "mode": mode,
                        "reactions": {
                            rid: {"before": {"lb": base[rid][0], "ub": base[rid][1], "kind": base[rid][2]},
                                  "after":  {"lb": after[rid][0], "ub": after[rid][1], "kind": after[rid][2]}}
                            for rid in after.keys()
                        }
                    })

                prog.advance(task)

    # write outputs
    stem = f"constraints__{mode}"
    tsv_path = outdir / f"{stem}.tsv"
    json_path = outdir / f"{stem}__reactions.json"

    # TSV summary
    if summary_rows:
        cols = ["spot_id", "microbe", "mode", "changed_ex", "changed_internal", "warnings"]
        with open(tsv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(summary_rows)
    else:
        # still write header
        with open(tsv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["spot_id","microbe","mode","changed_ex","changed_internal","warnings"])
            w.writeheader()

    # JSON full reaction bounds
    json_path.write_text(json.dumps(reactions_payload, indent=2))

    typer.secho(f"✅ Wrote {tsv_path.name} and {json_path.name} to {outdir.resolve()}", fg=typer.colors.GREEN)
