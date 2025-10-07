# microscape/cli/constrain.py
from __future__ import annotations
import json, csv, typer
from pathlib import Path
from datetime import datetime, timezone
from rich.progress import Progress
from typing import Dict, List, Tuple
from math import isfinite

from ..io.system_loader import load_system, iter_spot_files_for_env, read_spot_yaml, read_microbe_yaml
from ..io.metabolism_rules import load_rules, MetabolismRules, expr_to_exchange_bounds
from ..runner.constraints import constrain_one

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command("constrain")
def constrain_cmd(
    system_yml: Path = typer.Argument(..., help="Path to system.yml"),
    mode: str = typer.Option("combined", "--mode", case_sensitive=False,
                             help="Constraint source: environmental | transcriptional | combined"),
    outdir: Path = typer.Option("outputs/constraints", help="Output directory"),
    write_json: bool = typer.Option(True, "--json/--no-json", help="Write per-reaction JSON"),
    write_csv: bool = typer.Option(True, "--csv/--no-csv", help="Write summary TSV"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Extra logs"),
):
    """
    Build constraints for each spot×microbe and write:
    - summary TSV (one row per spot×microbe)
    - detailed per-reaction JSON (consumable by metabolism)
    """
    sys_info = load_system(system_yml)
    env_files = sys_info["environment_files"]
    outdir.mkdir(parents=True, exist_ok=True)

    # Resolve environmental mapping rules (your existing metabolism.yml use)
    metab_cfg_path = sys_info.get("metabolism_cfg")
    rules = load_rules(metab_cfg_path) if metab_cfg_path else None

    summary_rows: List[Dict] = []
    detail_nested: Dict = {
        "system": str(system_yml.resolve()),
        "mode": mode.lower(),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "solver": getattr(rules, "solver", None) if rules else None,
        "spots": {}
    }

    mode_l = mode.lower()
    if mode_l not in ("environmental","transcriptional","combined"):
        raise typer.BadParameter("mode must be one of: environmental | transcriptional | combined")

    def _rxn_type(r) -> str:
        # Fast path by naming convention
        if r.id.startswith("EX_"):
            return "exchange"

        mets = list(r.metabolites.keys())
        if len(mets) == 1:
            # Single-metabolite stoichiometry is typical for exchanges
            return "exchange"

        comps = {m.compartment for m in mets if hasattr(m, "compartment")}
        if "e" in comps and any(c != "e" for c in comps):
            return "transport"
        if comps == {"e"}:
            # Edge case: multi-stoich in extracellular only → treat as exchange
            return "exchange"
        return "internal"
    
    def get_base_bounds(model):
        out = {}
        for r in model.reactions:
            lb0 = float(r.lower_bound)
            ub0 = float(r.upper_bound)
            rtype = _rxn_type(r)
            out[r.id] = (lb0, ub0, rtype)
        return out

    def build_env_bounds(
        base_bounds: Dict[str, Tuple[float, float, str]],
        spot_mets: Dict[str, float],
    ) -> Dict[str, Tuple[Optional[float], Optional[float]]]:
        """
        Map extracellular metabolite concentrations (mM) to EX_ bounds using metabolism.yml rules.
        """
        env_bounds: Dict[str, Tuple[Optional[float], Optional[float]]] = {}
        if not rules:
            return env_bounds
        # rules.metabolite_map: { "<met id in spots>": "<EX_rxn_id>" }
        for met_id, ex_id in (rules.metabolite_map or {}).items():
            if ex_id not in base_bounds:
                continue
            # Only constrain exchanges
            lb0, ub0, rtype = base_bounds[ex_id]
            if rtype != "exchange":
                continue
            conc = float(spot_mets.get(met_id, 0.0))
            lb_env, ub_env = rules.uptake_to_bounds(conc)
            # Ensure lower ≤ upper
            if lb_env is not None and ub_env is not None and lb_env > ub_env:
                lb_env = ub_env
            env_bounds[ex_id] = (lb_env, ub_env)
        return env_bounds

    def build_tx_bounds(
        model,
        base_bounds: Dict[str, Tuple[float, float, str]],
        gene_tpm: Dict[str, float],
    ) -> Dict[str, Tuple[Optional[float], Optional[float]]]:
        """
        Convert transcript levels (TPM) → reaction "activity" via GPRs (AND=min, OR=max),
        then shrink INTERNAL reaction bounds toward zero proportionally to activity a∈[0,1].
          - irreversible forward (lb≥0): ub := ub*a
          - irreversible backward (ub≤0): lb := lb*a
          - reversible (lb<0<ub): lb := lb*a, ub := ub*a
        Exchanges are not touched here; they are governed by environmental rules.
        """
        tx_bounds: Dict[str, Tuple[Optional[float], Optional[float]]] = {}
        if not gene_tpm:
            return tx_bounds

        # Normalize TPMs per-spot to [0,1] for robustness without extra config
        max_tpm = max((float(v) for v in gene_tpm.values() if v is not None), default=0.0)
        gene_act = {g: (float(v) / max_tpm if max_tpm > 0 else 0.0) for g, v in gene_tpm.items() if v is not None}

        # Small recursive parser for rxn.gene_reaction_rule with AND/OR
        import re

        def gpr_activity(rule: str) -> float:
            rule = rule or ""
            if not rule.strip():
                return 1.0  # no GPR → don't constrain
            toks = re.findall(r"[A-Za-z0-9_]+|\(|\)|and|or", rule, flags=re.IGNORECASE)
            pos = 0

            def parse_expr():
                nonlocal pos
                vals = [parse_term()]
                while pos < len(toks) and toks[pos].lower() == "or":
                    pos += 1
                    vals.append(parse_term())
                return max(vals) if len(vals) > 1 else vals[0]

            def parse_term():
                nonlocal pos
                vals = [parse_factor()]
                while pos < len(toks) and toks[pos].lower() == "and":
                    pos += 1
                    vals.append(parse_factor())
                return min(vals) if len(vals) > 1 else vals[0]

            def parse_factor():
                nonlocal pos
                tok = toks[pos]
                pos += 1
                if tok == "(":
                    v = parse_expr()
                    assert toks[pos] == ")"
                    pos += 1
                    return v
                # gene symbol
                return float(gene_act.get(tok, 0.0))

            try:
                pos = 0
                v = parse_expr()
                return max(0.0, min(1.0, float(v)))
            except Exception:
                # If the rule can't be parsed, fail open (no extra constraint)
                return 1.0

        for rxn in model.reactions:
            rid = rxn.id
            lb0, ub0, rtype = base_bounds[rid]
            if rtype == "exchange":
                continue  # only internals here
            a = gpr_activity(rxn.gene_reaction_rule or "")
            if a >= 0.9999:
                continue  # no change
            # irreversible forward
            if lb0 >= 0.0 and ub0 > 0.0:
                tx_bounds[rid] = (lb0, ub0 * a)
            # irreversible backward
            elif ub0 <= 0.0 and lb0 < 0.0:
                tx_bounds[rid] = (lb0 * a, ub0)
            # reversible
            else:
                tx_bounds[rid] = (lb0 * a, ub0 * a)

        return tx_bounds

    with Progress() as prog:
        task = prog.add_task("[cyan]Constraining…", total=len(env_files))
        for env_file in env_files:
            for _, spot_path in iter_spot_files_for_env(env_file, sys_info["paths"]):
                spot = read_spot_yaml(spot_path)
                sid = spot.get("name") or spot.get("id") or spot_path.stem

                # Spot metabolite concentrations (for env constraints)
                met_vals = ((spot.get("measurements") or {}).get("metabolites") or {}).get("values") or {}

                # For each microbe we actually observe in this spot (or all in registry; pick your policy)
                microbes_vals = ((spot.get("measurements") or {}).get("microbes") or {}).get("values") or {}
                for mid in microbes_vals.keys():
                    myml = read_microbe_yaml(mid, sys_info)  # resolves via system paths
                    if not myml:
                        if verbose: typer.echo(f"WARN: Microbe YAML not found for {mid}")
                        continue
                    model_path = (Path(myml["__file__"]).parent / myml["microbe"]["model"]["path"]).resolve()

                    # Load model (cobra already installed in your env)
                    import cobra
                    try:
                        model = cobra.io.read_sbml_model(str(model_path))
                    except Exception as e:
                        if verbose:
                            typer.echo(f"ERROR loading SBML for {mid}: {e}")
                        continue

                    base_bounds = get_base_bounds(model)

                    # Environmental bounds from concentrations
                    env_bounds: Dict[str, Tuple[Optional[float], Optional[float]]] = {}
                    if mode_l in ("environmental", "combined") and rules:
                        env_bounds = build_env_bounds(base_bounds, met_vals)

                    # Transcriptional bounds from per-spot TPMs
                    tx_bounds: Dict[str, Tuple[Optional[float], Optional[float]]] = {}
                    if mode_l in ("transcriptional", "combined"):
                        tx_for_microbe = (
                            ((spot.get("measurements") or {}).get("transcripts") or {}).get("values") or {}
                        ).get(mid) or {}
                        tx_bounds = build_tx_bounds(model, base_bounds, tx_for_microbe)

                    # Merge to per-reaction decisions
                    summary_row, reaction_rows = constrain_one(
                        spot_id=sid,
                        microbe_id=mid,
                        base_bounds=base_bounds,
                        env_bounds=env_bounds,
                        tx_bounds=tx_bounds,
                    )
                    summary_row["mode"] = mode_l
                    summary_rows.append(summary_row)

                    # JSON assembly
                    dspot = detail_nested["spots"][env_id]["spots"].setdefault(sid, {"microbes": {}})
                    dmic = dspot["microbes"].setdefault(mid, {"reactions": {}, "summary": {}})
                    for rr in reaction_rows:
                        rid = rr["react_id"]
                        dmic["reactions"][rid] = {
                            "rtype": rr["rtype"],
                            "lb0": rr["lb0"],
                            "ub0": rr["ub0"],
                            "lb_env": rr["lb_env"],
                            "ub_env": rr["ub_env"],
                            "lb_tx": rr["lb_tx"],
                            "ub_tx": rr["ub_tx"],
                            "lb_final": rr["lb_final"],
                            "ub_final": rr["ub_final"],
                            "changed": rr["changed"],
                            "notes": rr["notes"],
                        }
                    dmic["summary"] = {
                        "changed_ex": summary_row["changed_ex"],
                        "changed_internal": summary_row["changed_internal"],
                        "warnings": summary_row["warnings"].split(";") if summary_row["warnings"] else [],
                    }
                prog.advance(task)

    # Write outputs
    stem = f"constraints__{mode_l}"
    if write_csv:
        tsv = outdir / f"{stem}.tsv"
        with tsv.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["spot_id","microbe","mode","changed_ex","changed_internal","warnings"])
            for r in summary_rows:
                w.writerow([r["spot_id"], r["microbe"], r["mode"], r["changed_ex"], r["changed_internal"], r["warnings"]])
        typer.echo(f"Summary TSV : {tsv}")

    if write_json:
        jpath = outdir / f"{stem}__reactions.json"
        jpath.write_text(json.dumps(detail_nested, indent=2))
        typer.echo(f"Detail JSON: {jpath}")

    typer.secho("✅ Constraints built.", fg=typer.colors.GREEN)
