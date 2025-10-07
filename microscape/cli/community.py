# microscape/cli/community.py
from __future__ import annotations
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional
import json
import csv
import math
import typer
from collections import defaultdict

app = typer.Typer(add_completion=False, no_args_is_help=True)

def _safe_float(x) -> float:
    try:
        return float(x)
    except Exception:
        return 0.0

def _split_exchange_fluxes(fluxes: Dict[str, Any]) -> Tuple[Dict[str, float], Dict[str, float]]:
    """
    Return (uptake, secretion) magnitude dicts keyed by exchange reaction id.
    COBRA sign convention: uptake is negative; secretion is positive.
    """
    uptake: Dict[str, float] = {}
    secretion: Dict[str, float] = {}
    for rid, v in (fluxes or {}).items():
        f = _safe_float(v)
        if f < 0:
            uptake[rid] = -f
        elif f > 0:
            secretion[rid] = f
    return uptake, secretion

def _schoener_D(uptake_i: Dict[str, float], uptake_j: Dict[str, float]) -> float:
    # proportions over union of resources
    keys = set(k for k, v in uptake_i.items() if v > 0) | set(k for k, v in uptake_j.items() if v > 0)
    if not keys:
        return 0.0
    sum_i = sum(uptake_i.get(k, 0.0) for k in keys)
    sum_j = sum(uptake_j.get(k, 0.0) for k in keys)
    if sum_i <= 0 or sum_j <= 0:
        return 0.0
    D = 1.0 - 0.5 * sum(abs(uptake_i.get(k, 0.0) / sum_i - uptake_j.get(k, 0.0) / sum_j) for k in keys)
    # clamp numeric noise
    return max(0.0, min(1.0, D))

def _jaccard_nonzero(uptake_i: Dict[str, float], uptake_j: Dict[str, float]) -> float:
    A = {k for k, v in uptake_i.items() if v > 0}
    B = {k for k, v in uptake_j.items() if v > 0}
    if not A and not B:
        return 0.0
    return len(A & B) / float(len(A | B)) if (A or B) else 0.0

def _cross_feeding(secretion_i: Dict[str, float], uptake_j: Dict[str, float]) -> float:
    # same exchange reaction id is used for secreted resource consumed by j
    keys = set(secretion_i.keys()) & set(uptake_j.keys())
    return sum(min(secretion_i[k], uptake_j[k]) for k in keys)

def _unique_resource_fraction(uptake_i: Dict[str, float], others: List[Dict[str, float]]) -> float:
    total_i = sum(uptake_i.values())
    if total_i <= 0:
        return 0.0
    others_union = set()
    for up in others:
        others_union |= {k for k, v in up.items() if v > 0}
    unique_keys = [k for k, v in uptake_i.items() if v > 0 and k not in others_union]
    unique_uptake = sum(uptake_i[k] for k in unique_keys)
    return unique_uptake / total_i

def _redundancy(uptakes: List[Dict[str, float]]) -> float:
    """
    Resource redundancy across the community: average fraction of microbes that use each resource.
    Values in [0,1]. Higher = more shared resource use (more competitive pressure).
    """
    if not uptakes:
        return 0.0
    N = len(uptakes)
    all_keys = set().union(*[{k for k, v in up.items() if v > 0} for up in uptakes])
    if not all_keys:
        return 0.0
    frac_sum = 0.0
    for k in all_keys:
        users = sum(1 for up in uptakes if up.get(k, 0.0) > 0)
        frac_sum += users / N
    return frac_sum / len(all_keys)

def _connectance_cf(cf_mat: Dict[str, Dict[str, float]]) -> float:
    """
    Directed connectance of cross-feeding graph: fraction of possible i->j (i!=j) with CF_ij>0.
    """
    microbes = list(cf_mat.keys())
    if len(microbes) <= 1:
        return 0.0
    M = len(microbes)
    possible = M * (M - 1)
    edges = 0
    for i in microbes:
        for j in microbes:
            if i == j:
                continue
            if cf_mat[i].get(j, 0.0) > 0:
                edges += 1
    return edges / possible

def _mean(xs: List[float]) -> float:
    return sum(xs) / len(xs) if xs else 0.0

def _median(xs: List[float]) -> float:
    if not xs:
        return 0.0
    ys = sorted(xs)
    n = len(ys)
    mid = n // 2
    return (ys[mid] if n % 2 else 0.5 * (ys[mid - 1] + ys[mid]))

def _load_metabolism_json(path: Path) -> Dict[str, Any]:
    data = json.loads(path.read_text())
    # expected shape: { "<spot_id>": { "microbes": { "<mid>": {"fluxes": {...}, "objective": float, ...}, ... } } }
    return data

@app.command("community")
def community_cmd(
    metabolism_json: Path = typer.Argument(..., help="Path to a named metabolism JSON (e.g., metabolism_constrained_transcriptional.json)"),
    outdir: Path = typer.Option("outputs/community", help="Output directory for community metrics"),
    write_graph: bool = typer.Option(False, "--graph/--no-graph", help="Also write community_graph.json with bipartite & cross-feeding edges"),
):
    """
    Compute community-level metrics (competition, complementarity, cross-feeding)
    from named metabolism outputs.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    per_spot = _load_metabolism_json(metabolism_json)

    rows_spot: List[Dict[str, Any]] = []
    rows_pair: List[Dict[str, Any]] = []
    rows_microbe: List[Dict[str, Any]] = []
    graph = {"spots": {}} if write_graph else None

    for spot_id, srec in per_spot.items():
        microbes: Dict[str, Any] = (srec or {}).get("microbes", {}) or {}
        if not microbes:
            continue

        # Build uptake/secretion dicts per microbe
        uptake_by_mid: Dict[str, Dict[str, float]] = {}
        secre_by_mid: Dict[str, Dict[str, float]] = {}
        obj_by_mid: Dict[str, float] = {}

        for mid, mrec in microbes.items():
            fluxes = (mrec or {}).get("fluxes", {}) or {}
            uptake, secre = _split_exchange_fluxes(fluxes)
            uptake_by_mid[mid] = uptake
            secre_by_mid[mid] = secre
            obj_by_mid[mid] = _safe_float((mrec or {}).get("objective", 0))

        mids = sorted(uptake_by_mid.keys())
        if len(mids) == 0:
            continue

        # Per-microbe summary (unique resource use, totals)
        for mid in mids:
            others = [uptake_by_mid[o] for o in mids if o != mid]
            urf = _unique_resource_fraction(uptake_by_mid[mid], others)
            rows_microbe.append({
                "spot_id": spot_id,
                "microbe": mid,
                "unique_resource_fraction": urf,
                "total_uptake": sum(uptake_by_mid[mid].values()),
                "total_secretion": sum(secre_by_mid[mid].values()),
                "objective": obj_by_mid.get(mid, 0.0),
            })

        # Pairwise competition/complementarity + cross-feeding
        comp_vals: List[float] = []
        compl_vals: List[float] = []
        cf_mat: Dict[str, Dict[str, float]] = {i: {} for i in mids}

        for i in mids:
            for j in mids:
                if i == j:
                    continue
                D = _schoener_D(uptake_by_mid[i], uptake_by_mid[j])   # competition similarity
                C = 1.0 - D                                           # complementarity
                CF_ij = _cross_feeding(secre_by_mid[i], uptake_by_mid[j])
                cf_mat[i][j] = CF_ij

                if i < j:
                    comp_vals.append(D)
                    compl_vals.append(C)
                rows_pair.append({
                    "spot_id": spot_id,
                    "microbe_i": i,
                    "microbe_j": j,
                    "competition_schoenerD": D,
                    "complementarity": C,
                    "crossfeeding_ij": CF_ij,
                    # also add reverse if present to help spotting asymmetry
                    "crossfeeding_ji": None,  # filled in a second pass
                    "net_crossfeeding_i_minus_j": None,
                    "uptake_jaccard": _jaccard_nonzero(uptake_by_mid[i], uptake_by_mid[j]),
                })

        # Fill reverse CF and net on the pair table
        # Index for faster lookup
        pair_index = {}
        for r in rows_pair:
            if r["spot_id"] == spot_id:
                pair_index[(r["microbe_i"], r["microbe_j"])] = r
        for i in mids:
            for j in mids:
                if i == j:
                    continue
                CF_ij = cf_mat[i].get(j, 0.0)
                CF_ji = cf_mat[j].get(i, 0.0)
                rec = pair_index.get((i, j))
                if rec is not None:
                    rec["crossfeeding_ji"] = CF_ji
                    rec["net_crossfeeding_i_minus_j"] = CF_ij - CF_ji

        # Spot-level rollups
        rows_spot.append({
            "spot_id": spot_id,
            "n_microbes": len(mids),
            "mean_competition": _mean(comp_vals),
            "median_competition": _median(comp_vals),
            "mean_complementarity": _mean(compl_vals),
            "median_complementarity": _median(compl_vals),
            "resource_redundancy": _redundancy([uptake_by_mid[m] for m in mids]),
            "crossfeeding_connectance": _connectance_cf(cf_mat),
            "total_uptake_sum": sum(sum(up.values()) for up in uptake_by_mid.values()),
            "total_secretion_sum": sum(sum(sc.values()) for sc in secre_by_mid.values()),
            "objective_sum": sum(obj_by_mid.values()),
        })

        # Optional graph
        if write_graph:
            # bipartite microbe <-> exchange; and derived microbe->microbe edges weighted by CF
            rnodes = set()
            for mid in mids:
                rnodes |= set(uptake_by_mid[mid].keys()) | set(secre_by_mid[mid].keys())
            g_spot = {"microbes": mids, "reactions": sorted(rnodes), "edges": {"uptake": [], "secretion": [], "crossfeeding": []}}
            for mid in mids:
                for rid, mag in uptake_by_mid[mid].items():
                    if mag > 0:
                        g_spot["edges"]["uptake"].append([mid, rid, mag])
                for rid, mag in secre_by_mid[mid].items():
                    if mag > 0:
                        g_spot["edges"]["secretion"].append([mid, rid, mag])
            for i in mids:
                for j in mids:
                    if i != j and cf_mat[i].get(j, 0.0) > 0:
                        g_spot["edges"]["crossfeeding"].append([i, j, cf_mat[i][j]])
            graph["spots"][spot_id] = g_spot

    # Write outputs
    outdir.mkdir(parents=True, exist_ok=True)

    p_spot = outdir / "community_metrics_spot.csv"
    with p_spot.open("w", newline="") as fh:
        cols = ["spot_id","n_microbes","mean_competition","median_competition",
                "mean_complementarity","median_complementarity",
                "resource_redundancy","crossfeeding_connectance",
                "total_uptake_sum","total_secretion_sum","objective_sum"]
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows_spot:
            w.writerow(r)

    p_pair = outdir / "community_pairwise.csv"
    with p_pair.open("w", newline="") as fh:
        cols = ["spot_id","microbe_i","microbe_j",
                "competition_schoenerD","complementarity",
                "uptake_jaccard","crossfeeding_ij","crossfeeding_ji",
                "net_crossfeeding_i_minus_j"]
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows_pair:
            w.writerow(r)

    p_micro = outdir / "community_microbe.csv"
    with p_micro.open("w", newline="") as fh:
        cols = ["spot_id","microbe","unique_resource_fraction","total_uptake","total_secretion","objective"]
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows_microbe:
            w.writerow(r)

    if write_graph:
        (outdir / "community_graph.json").write_text(json.dumps(graph, indent=2))

    typer.secho("✅ Community metrics written.", fg=typer.colors.GREEN)
    typer.echo(f"  Spot summary  : {p_spot}")
    typer.echo(f"  Pairwise      : {p_pair}")
    typer.echo(f"  Microbe-level : {p_micro}")
    if write_graph:
        typer.echo(f"  Graph JSON    : {(outdir / 'community_graph.json')}")
