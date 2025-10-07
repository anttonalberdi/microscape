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


def _select_transcripts(spot: dict, mid: str, sid: str) -> Dict[str, float]:
    """
    Robustly select a gene->TPM dict for this spot×microbe.
    Supports:
      1) per-microbe dict: values[mid] = {gene: tpm, ...}
      2) shared per-spot dict: values = {gene: tpm, ...}
      3) nested dicts under keys like 'all', spot id/name, etc.
    Returns {} if nothing sensible is found.
    """
    tx_root = ((spot.get("measurements") or {}).get("transcripts") or {}).get("values") or {}
    if not isinstance(tx_root, dict):
        return {}

    # 1) direct per-microbe
    v = tx_root.get(mid)
    if isinstance(v, dict) and v:
        return v

    # 2) flat dict (shared for everyone)
    # heuristic: all values are scalar-like (numbers/str)
    if tx_root and all(not isinstance(x, dict) for x in tx_root.values()):
        return tx_root  # assume {gene: TPM}

    # 3) look for common keys then fallback to first dict value
    for key in (sid, spot.get("name"), spot.get("id"), "all", "ALL", "default", "shared"):
        if key and isinstance(tx_root.get(key), dict) and tx_root.get(key):
            return tx_root[key]

    for _k, _v in tx_root.items():
        if isinstance(_v, dict) and _v:
            return _v

    return {}


def _build_tx_bounds(
    model,
    base_bounds: Dict[str, Tuple[float, float, str]],
    gene_tpm_in: Dict[str, float],
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
    if not gene_tpm_in:
        return tx_bounds

    # case-insensitive gene matching (GPR tokens vs TPM keys)
    # keep numeric conversion robust
    gene_tpm = {}
    for g, v in gene_tpm_in.items():
        try:
            gene_tpm[str(g).lower()] = float(v)
        except Exception:
            continue

    if not gene_tpm:
        return tx_bounds

    max_tpm = max(gene_tpm.values()) if gene_tpm else 0.0
    gene_act = {g: (val / max_tpm if max_tpm > 0 else 0.0) for g, val in gene_tpm.items()}

    import re

    def gpr_activity(rule: str) -> float:
        rule = rule or ""
        if not rule.strip():
            return 1.0
        # allow underscores, dots, and hyphens in gene IDs
        toks = re.findall(r"[A-Za-z0-9_.-]+|\(|\)|and|or", rule, flags=re.IGNORECASE)
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
            # gene symbol; look up case-insensitively
            return float(gene_act.get(tok.lower(), 0.0))

        try:
            pos = 0
            v = parse_expr()
            return max(0.0, min(1.0, float(v)))
        except Exception:
            return 1.0  # fail open (do not over-constrain)

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
    write_json: bool = typer.Option(True, "--json/--no-json", help="Write JSON outputs"),
    write_csv: bool = typer.Option(True, "--csv/--no-csv", help="Write summary CSV"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Extra logs"),
):
    """
    Build constraints for each spot×microbe and write:
      - summary CSV (one row per spot×microbe)
      - detailed JSON (debug)
      - compact constraints JSON (for `microscape metabolism --constraints`)
    """
    sys_info = load_system(system_yml)
    env_files = sys_info["environment_files"]
    outdir.mkdir(parents=True, exist_ok=True)

    metab_cfg_path = sys_info.get("metabolism_cfg")
    rules: Optional[MetabolismRules] = load_rules(metab_cfg_path) if metab_cfg_path else None

    summary_rows: List[Dict] = []

    # detailed/debug JSON (rich info)
    debug_nested: Dict = {
        "system": str(system_yml.resolve()),
        "mode": mode.lower(),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "spots": {},
    }

    # compact JSON strictly for metabolism (final bounds only; changed reactions)
    metab_nested: Dict = {
        "version": 1,
        "system": str(system_yml.resolve()),
        "mode": mode.lower(),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "spots": {}  # spots[spot_id].microbes[microbe_id].reactions[rxn_id] = {"lb": .., "ub": ..}
    }

    mode_l = mode.lower()
    if mode_l not in ("environmental", "transcriptional", "combined"):
        raise typer.BadParameter("mode must be one of: environmental | transcriptional | combined")

    with Progress() as prog:
        task = prog.add_task("[cyan]Constraining…", total=len(env_files))
        for env_file in env_files:
            # Spot discovery anchored to system.yml's paths.spots_dir
            spot_pairs = _iter_spots_for_env_via_system(env_file, sys_info)
            env_doc = yaml.safe_load(Path(env_file).read_text()) or {}
            env_id = (env_doc.get("environment") or {}).get("id") or Path(env_file).stem

            debug_nested["spots"].setdefault(env_id, {"spots": {}})

            for sid, spot_path in spot_pairs:
                spot = read_spot_yaml(spot_path) or {}
                sid = spot.get("name") or spot.get("id") or sid

                met_vals = ((spot.get("measurements") or {}).get("metabolites") or {}).get("values") or {}
                microbes_vals = ((spot.get("measurements") or {}).get("microbes") or {}).get("values") or {}

                # set up metab JSON node
                metab_spot = metab_nested["spots"].setdefault(sid, {"microbes": {}})

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
                        # ✅ robust transcript selection
                        gene_tpm = _select_transcripts(spot, mid, sid)
                        tx_bounds = _build_tx_bounds(model, base_bounds, gene_tpm)
                        if verbose and not gene_tpm:
                            typer.echo(f"INFO: no transcripts found for {sid} × {mid}; internal bounds unchanged")

                    summary_row, reaction_rows = constrain_one(
                        spot_id=sid,
                        microbe_id=mid,
                        base_bounds=base_bounds,
                        env_bounds=env_bounds,
                        tx_bounds=tx_bounds,
                    )
                    summary_row["mode"] = mode_l
                    summary_rows.append(summary_row)

                    # ----- detailed debug JSON -----
                    dspot = debug_nested["spots"][env_id]["spots"].setdefault(sid, {"microbes": {}})
                    dmic = dspot["microbes"].setdefault(mid, {"reactions": {}, "summary": {}})
                    # also prepare metab JSON node
                    mnode = metab_spot["microbes"].setdefault(mid, {"reactions": {}})

                    for rr in reaction_rows:
                        rid = rr.get("react_id") or rr.get("reaction") or rr.get("id")
                        if rid is None:
                            continue
                        rr_type = rr.get("rtype") or rr.get("type") or ("exchange" if str(rid).startswith("EX_") else "internal")

                        # rich record for humans
                        dmic["reactions"][rid] = {
                            "type": rr_type,
                            "lb0": rr.get("lb0"),
                            "ub0": rr.get("ub0"),
                            "lb_env": rr.get("lb_env"),
                            "ub_env": rr.get("ub_env"),
                            "lb_tx": rr.get("lb_tx"),
                            "ub_tx": rr.get("ub_tx"),
                            "lb_final": rr.get("lb_final"),
                            "ub_final": rr.get("ub_final"),
                            "changed": rr.get("changed", False),
                            "notes": rr.get("notes", ""),
                        }
                        # compact record for metabolism (only changed reactions)
                        if rr.get("changed", False):
                            mnode["reactions"][rid] = {
                                "lb": rr.get("lb_final"),
                                "ub": rr.get("ub_final"),
                            }

                    dmic["summary"] = {
                        "changed_ex": summary_row["changed_ex"],
                        "changed_internal": summary_row["changed_internal"],
                        "warnings": summary_row["warnings"].split(";") if summary_row["warnings"] else [],
                    }

            prog.advance(task)

    # ---------- write outputs ----------
    stem = f"constraints__{mode_l}"

    if write_csv:
        csv_path = outdir / f"{stem}.csv"
        with csv_path.open("w", newline="") as fh:
            w = csv.writer(fh)  # default delimiter = ','
            w.writerow(["spot_id","microbe","mode","changed_ex","changed_internal","warnings"])
            for r in summary_rows:
                w.writerow([r["spot_id"], r["microbe"], r["mode"], r["changed_ex"], r["changed_internal"], r["warnings"]])
        typer.echo(f"Summary CSV : {csv_path}")

    if write_json:
        # Detailed debug JSON (keeps rich info you already had)
        j_debug = outdir / f"{stem}__debug.json"
        j_debug.write_text(json.dumps(debug_nested, indent=2))
        typer.echo(f"Debug JSON  : {j_debug}")

        # Compact constraints JSON for metabolism
        j_metab = outdir / f"{stem}__reactions.json"
        j_metab.write_text(json.dumps(metab_nested, indent=2))
        typer.echo(f"Metab JSON  : {j_metab}")

    typer.secho("✅ Constraints built.", fg=typer.colors.GREEN)
