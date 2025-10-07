# microscape/cli/constrain.py
from __future__ import annotations
import json, csv, typer, yaml
from pathlib import Path
from datetime import datetime, timezone
from rich.progress import Progress
from typing import Dict, List, Tuple, Optional

from ..io.system_loader import (
    load_system,
    read_spot_yaml,
    read_microbe_yaml,
)
from ..io.metabolism_rules import load_rules, MetabolismRules
from ..runner.constraints import constrain_one

app = typer.Typer(add_completion=False, no_args_is_help=True)


# ------- helpers -------

def _rxn_type(r) -> str:
    if r.id.startswith("EX_"):
        return "exchange"
    mets = list(r.metabolites.keys())
    if len(mets) == 1:
        return "exchange"
    comps = {getattr(m, "compartment", None) for m in mets}
    if "e" in comps and any(c != "e" for c in comps):
        return "transport"
    if comps == {"e"}:
        return "exchange"
    return "internal"


def _get_base_bounds(model):
    out = {}
    for r in model.reactions:
        lb0 = float(r.lower_bound)
        ub0 = float(r.upper_bound)
        out[r.id] = (lb0, ub0, "exchange" if _rxn_type(r) == "exchange" else "internal")
    return out


def _iter_spots_for_env_via_system(env_file: Path, sys_info: Dict) -> List[Tuple[str, Path]]:
    """
    Return [(spot_id, spot_path), ...] for the given environment file,
    by scanning the spots_dir configured in system.yml and selecting those
    whose spot.env_id matches this environment's id.
    This avoids any accidental search under the environments directory.
    """
    env_doc = yaml.safe_load(Path(env_file).read_text()) or {}
    env = env_doc.get("environment") or {}
    env_id = env.get("id") or Path(env_file).stem

    sys_paths = sys_info.get("paths") or {}
    spots_dir = Path(sys_paths.get("spots_dir") or "")
    if not spots_dir.is_absolute():
        root = Path(sys_info.get("root") or Path(sys_info.get("system_path", ".")).parent)
        spots_dir = (root / spots_dir).resolve()

    out: List[Tuple[str, Path]] = []
    candidates = sorted(list(spots_dir.glob("*.yml")) + list(spots_dir.glob("*.yaml")))
    for p in candidates:
        try:
            raw = yaml.safe_load(p.read_text()) or {}
        except Exception:
            continue
        s = (raw.get("spot") or raw)
        sid = s.get("id") or s.get("name") or p.stem
        env_in_spot = s.get("env_id") or s.get("environment_id")
        if env_in_spot == env_id:
            out.append((sid, p))
    return out


def _build_env_bounds(
    rules: Optional[MetabolismRules],
    base_bounds: Dict[str, Tuple[float, float, str]],
    spot_mets: Dict[str, float],
) -> Dict[str, Tuple[Optional[float], Optional[float]]]:
    """Map extracellular metabolite concentrations (mM) to EX_ bounds using metabolism.yml rules."""
    env_bounds: Dict[str, Tuple[Optional[float], Optional[float]]] = {}
    if not rules:
        return env_bounds
    for met_id, ex_id in (rules.metabolite_map or {}).items():
        if ex_id not in base_bounds:
            continue
        lb0, ub0, rtype = base_bounds[ex_id]
        if rtype != "exchange":
            continue
        conc = float(spot_mets.get(met_id, 0.0) or 0.0)
        lb_env, ub_env = rules.uptake_to_bounds(conc)
        if lb_env is not None and ub_env is not None and lb_env > ub_env:
            lb_env = ub_env
        env_bounds[ex_id] = (lb_env, ub_env)
    return env_bounds


def _build_tx_bounds(
    model,
    base_bounds: Dict[str, Tuple[float, float, str]],
    gene_tpm: Dict[str, float],
) -> Dict[str, Tuple[Optional[float], Optional[float]]]:
    """
    TPM → reaction activity via GPRs (AND=min, OR=max).
    Tighten INTERNAL reaction bounds toward zero:
      - irreversible forward (lb≥0): ub := ub*a
      - irreversible backward (ub≤0): lb := lb*a
      - reversible (lb<0<ub): lb := lb*a, ub := ub*a
    Exchanges are governed by environmental rules.
    """
    tx_bounds: Dict[str, Tuple[Optional[float], Optional[float]]] = {}
    if not gene_tpm:
        return tx_bounds

    max_tpm = max((float(v) for v in gene_tpm.values() if v is not None), default=0.0)
    gene_act = {g: (float(v) / max_tpm if max_tpm > 0 else 0.0) for g, v in gene_tpm.items() if v is not None}

    import re

    def gpr_activity(rule: str) -> float:
        rule = rule or ""
        if not rule.strip():
            return 1.0
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
            return float(gene_act.get(tok, 0.0))

        try:
            pos = 0
            v = parse_expr()
            return max(0.0, min(1.0, float(v)))
        except Exception:
            return 1.0

    for rxn in model.reactions:
        rid = rxn.id
        lb0, ub0, rtype = base_bounds[rid]
        if rtype == "exchange":
            continue
        a = gpr_activity(rxn.gene_reaction_rule or "")
        if a >= 0.9999:
            continue
        if lb0 >= 0.0 and ub0 > 0.0:
            tx_bounds[rid] = (lb0, ub0 * a)
        elif ub0 <= 0.0 and lb0 < 0.0:
            tx_bounds[rid] = (lb0 * a, ub0)
        else:
            tx_bounds[rid] = (lb0 * a, ub0 * a)

    return tx_bounds


# ------- CLI -------

@app.command("constrain")
def constrain_cmd(
    system_yml: Path = typer.Argument(..., help="Path to system.yml"),
    mode: str = typer.Option(
        "combined",
        "--mode",
        case_sensitive=False,
        help="Constraint source: environmental | transcriptional | combined",
    ),
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

    metab_cfg_path = sys_info.get("metabolism_cfg")
    rules: Optional[MetabolismRules] = load_rules(metab_cfg_path) if metab_cfg_path else None

    summary_rows: List[Dict] = []
    detail_nested: Dict = {
        "system": str(system_yml.resolve()),
        "mode": mode.lower(),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "spots": {},
    }

    mode_l = mode.lower()
    if mode_l not in ("environmental", "transcriptional", "combined"):
        raise typer.BadParameter("mode must be one of: environmental | transcriptional | combined")

    with Progress() as prog:
        task = prog.add_task("[cyan]Constraining…", total=len(env_files))
        for env_file in env_files:
            # ✅ FIX: enumerate spots from system.yml's spots_dir, filtered by env_id
            spot_pairs = _iter_spots_for_env_via_system(env_file, sys_info)
            env_doc = yaml.safe_load(Path(env_file).read_text()) or {}
            env_id = (env_doc.get("environment") or {}).get("id") or Path(env_file).stem
            detail_nested["spots"].setdefault(env_id, {"spots": {}})

            for sid, spot_path in spot_pairs:
                spot = read_spot_yaml(spot_path) or {}
                sid = spot.get("name") or spot.get("id") or sid

                met_vals = ((spot.get("measurements") or {}).get("metabolites") or {}).get("values") or {}
                microbes_vals = ((spot.get("measurements") or {}).get("microbes") or {}).get("values") or {}

                for mid in microbes_vals.keys():
                    myml = read_microbe_yaml(mid, sys_info)
                    if not myml:
                        if verbose:
                            typer.echo(f"WARN: Microbe YAML not found for {mid}")
                        continue
                    model_path = (Path(myml["__file__"]).parent / myml["microbe"]["model"]["path"]).resolve()

                    import cobra
                    try:
                        model = cobra.io.read_sbml_model(str(model_path))
                    except Exception as e:
                        if verbose:
                            typer.echo(f"ERROR loading SBML for {mid}: {e}")
                        continue

                    base_bounds = _get_base_bounds(model)

                    env_bounds: Dict[str, Tuple[Optional[float], Optional[float]]] = {}
                    if mode_l in ("environmental", "combined") and rules:
                        env_bounds = _build_env_bounds(rules, base_bounds, met_vals)

                    tx_bounds: Dict[str, Tuple[Optional[float], Optional[float]]] = {}
                    if mode_l in ("transcriptional", "combined"):
                        tx_for_microbe = (
                            ((spot.get("measurements") or {}).get("transcripts") or {}).get("values") or {}
                        ).get(mid) or {}
                        tx_bounds = _build_tx_bounds(model, base_bounds, tx_for_microbe)

                    summary_row, reaction_rows = constrain_one(
                        spot_id=sid,
                        microbe_id=mid,
                        base_bounds=base_bounds,
                        env_bounds=env_bounds,
                        tx_bounds=tx_bounds,
                    )
                    summary_row["mode"] = mode_l
                    summary_rows.append(summary_row)

                    dspot = detail_nested["spots"][env_id]["spots"].setdefault(sid, {"microbes": {}})
                    dmic = dspot["microbes"].setdefault(mid, {"reactions": {}, "summary": {}})
                    for rr in reaction_rows:
                        rid = rr.pop("react_id")
                        rtype = rr.pop("type")
                        entry = {
                            "type": rtype,
                            "lb_env": rr["lb_env"],
                            "ub_env": rr["ub_env"],
                            "lb_tx": rr["lb_tx"],
                            "ub_tx": rr["ub_tx"],
                            "lb_final": rr["lb_final"],
                            "ub_final": rr["ub_final"],
                            "changed": rr["changed"],
                            "notes": rr["notes"],
                        }
                        dmic["reactions"][rid] = entry
                    dmic["summary"] = {
                        "changed_ex": summary_row["changed_ex"],
                        "changed_internal": summary_row["changed_internal"],
                        "warnings": summary_row["warnings"].split(";") if summary_row["warnings"] else [],
                    }

            prog.advance(task)

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
