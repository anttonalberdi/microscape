from __future__ import annotations
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
import json
import math
import numpy as np
import pandas as pd
import typer
import yaml

try:
    # Optional; if present we use CSR to write adjacencies
    from scipy.sparse import csr_matrix, save_npz
except Exception:
    csr_matrix = None
    save_npz = None

app = typer.Typer(add_completion=False, no_args_is_help=True)

# ───────────────────────────── helpers ─────────────────────────────

def _safe_float(x):
    try:
        return float(x)
    except Exception:
        return None

def _load_yaml(p: Path) -> dict:
    return yaml.safe_load(p.read_text()) or {}

def _iter_env_files(system_yml: Path) -> List[Path]:
    doc = _load_yaml(system_yml)
    sysn = (doc.get("system") or {})
    envs = sysn.get("environments") or []
    base = system_yml.parent
    return [(base / e).resolve() for e in envs]

def _iter_spots_for_env(env_file: Path, system_root: Path) -> List[Tuple[str, Path]]:
    env = _load_yaml(env_file)
    e = env.get("environment") or {}
    spots = e.get("spots") or []
    out: List[Tuple[str, Path]] = []
    for s in spots:
        if isinstance(s, str):
            sid = Path(s).stem
            # allow "spots/S0001.yml" or just "S0001.yml"
            sp = (system_root / s).resolve()
            if not sp.exists():
                sp = (system_root / "spots" / s).resolve()
        else:
            sid = s.get("id") or s.get("name") or Path(s.get("path")).stem
            sp = (system_root / s.get("path")).resolve()
        out.append((str(sid), sp))
    return out

def _load_spot(spot_file: Path) -> Dict[str, Any]:
    return _load_yaml(spot_file)

def _compute_neighbors_xy(coords: List[Tuple[float,float]], radius: Optional[float]=None, k: Optional[int]=None) -> List[List[int]]:
    """
    Brute-force neighbour lists for each point.
    If radius given: points within radius (including self).
    If k given: k nearest (including self).
    If neither: [i] only.
    """
    n = len(coords)
    if n == 0:
        return []
    xs = np.array([c[0] for c in coords], dtype=float)
    ys = np.array([c[1] for c in coords], dtype=float)
    valid = np.isfinite(xs) & np.isfinite(ys)
    neigh: List[List[int]] = [[] for _ in range(n)]
    if radius is None and (k is None or k <= 1):
        for i in range(n):
            neigh[i] = [i]
        return neigh
    for i in range(n):
        if not valid[i]:
            neigh[i] = []
            continue
        di: List[Tuple[float,int]] = []
        xi, yi = xs[i], ys[i]
        for j in range(n):
            if not valid[j]:
                continue
            dx, dy = xi - xs[j], yi - ys[j]
            d = math.hypot(dx, dy)
            if radius is not None:
                if d <= float(radius):
                    di.append((d, j))
            else:
                di.append((d, j))
        if k is not None:
            di = sorted(di, key=lambda t: t[0])[:max(1, int(k))]
        neigh[i] = [j for _, j in di]
    return neigh

def _adjacency_from_neigh(neigh: List[List[int]], symmetrize: bool=True, self_loops: bool=False):
    n = len(neigh)
    rows, cols, data = [], [], []
    for i, js in enumerate(neigh):
        for j in js:
            if (not self_loops) and (i == j):
                continue
            rows.append(i); cols.append(j); data.append(1)
    if csr_matrix is None:
        # Minimal dense fallback (only used if SciPy missing)
        import numpy as _np
        A = _np.zeros((n, n), dtype=int)
        for r, c in zip(rows, cols):
            A[r, c] = 1
        if symmetrize:
            A = ((A + A.T) > 0).astype(int)
        return A
    A = csr_matrix((data, (rows, cols)), shape=(n, n))
    if symmetrize:
        A = ((A + A.T) > 0).astype(int)
    return A

def _index_metabolism(meta: Dict[str, Any]) -> Dict[str, Any]:
    """Normalize metabolism JSON into spots -> microbes dictionary."""
    if not meta:
        return {}
    if "spots" in meta and isinstance(meta["spots"], dict):
        return meta["spots"]
    return meta  # assume already spot keyed

def _read_metabolism(metabolism_json: Optional[Path]) -> Dict[str, Any]:
    if not metabolism_json:
        return {}
    p = Path(metabolism_json)
    if not p.exists():
        raise typer.BadParameter(f"Metabolism JSON not found: {p}")
    try:
        return json.loads(p.read_text())
    except Exception as e:
        raise typer.BadParameter(f"Could not read metabolism JSON: {p} ({e})")

def _get_target_value(rec: Dict[str, Any], selector: str) -> Optional[float]:
    """
    selector:
      - 'objective'
      - 'flux:EX_Cxxxx_e'
    rec: microbe record from metabolism JSON
    """
    if selector == "objective":
        return _safe_float(rec.get("objective"))
    if selector.startswith("flux:"):
        rid = selector.split(":", 1)[1]
        fx = rec.get("fluxes") or rec.get("fluxes_all_exchanges") or {}
        return _safe_float(fx.get(rid))
    return None

def _combine_targets(values: List[Optional[float]], op: str) -> Optional[float]:
    vals = [v for v in values if v is not None and np.isfinite(v)]
    if not vals:
        return None
    op = (op or "sum").lower()
    if op in ("sum",):
        return float(np.sum(vals))
    if op in ("mean", "avg"):
        return float(np.mean(vals))
    if op in ("min",):
        return float(np.min(vals))
    if op in ("max",):
        return float(np.max(vals))
    if op in ("uptake_sum", "uptake"):
        # sum of negative flux magnitudes
        return float(np.sum([-min(0.0, v) for v in vals]))
    if op in ("secretion_sum", "secretion"):
        return float(np.sum([max(0.0, v) for v in vals]))
    # default to sum
    return float(np.sum(vals))

# ───────────────────────────── CLI ─────────────────────────────

@app.command("build")
def build_cmd(
    system_yml: Path = typer.Argument(..., help="Path to system.yml"),
    outdir: Path = typer.Option("outputs/model", help="Directory for outputs."),
    # Target assembly (unchanged capabilities)
    metabolism_json: Optional[Path] = typer.Option(None, "--metabolism-json", help="Metabolism JSON to pull targets from."),
    target: List[str] = typer.Option(None, "--target", help="Target selector(s): 'objective' or 'flux:EX_Cxxxx_e'. Repeatable."),
    target_op: Optional[str] = typer.Option(None, "--target-op", help="Combine multiple --target values: sum|mean|min|max|uptake_sum|secretion_sum."),
    # Spatial feature options (all optional / opt-in)
    include_coords: bool = typer.Option(False, help="Include x_um,y_um,z_um as features."),
    coords_poly: bool = typer.Option(False, help="Add x_um^2, y_um^2, x_um*y_um."),
    spatial_agg_radius: Optional[float] = typer.Option(None, help="Neighbourhood radius (µm) for mean features."),
    spatial_agg_k: Optional[int] = typer.Option(None, help="Neighbourhood k for kNN mean features."),
):
    """
    Assemble the modelling table with original target logic (multiple --target + --target-op)
    and optional spatial features (coords, polynomial terms, neighbourhood aggregates).
    """
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    spatial_dir = outdir / "spatial"
    spatial_dir.mkdir(exist_ok=True)

    # Load metabolism (for targets)
    metab = _read_metabolism(metabolism_json)
    metab_spots = _index_metabolism(metab)

    # Validate target usage (backward compatible: if user passed --target, we require it to build 'target' column)
    if target and not metabolism_json:
        raise typer.BadParameter("You passed --target but did not provide --metabolism-json.")
    if (not target) and metabolism_json:
        # allowed: the user might only want features table; we won't create a 'target' column
        pass

    system_root = Path(system_yml).parent.resolve()
    env_files = _iter_env_files(system_yml)

    rows: List[Dict[str, Any]] = []

    # Iterate environments
    for env_file in env_files:
        env_doc = _load_yaml(env_file)
        env_meta = (env_doc.get("environment") or {})
        env_id = str(env_meta.get("id") or Path(env_file).stem)

        # collect env spots (order matters for adjacency)
        spot_ids: List[str] = []
        spot_paths: List[Path] = []
        coords_xy: List[Tuple[float,float]] = []

        for sid, spath in _iter_spots_for_env(env_file, system_root):
            sobj = _load_spot(spath) or {}
            s = sobj.get("spot") or sobj
            spot_id = str(s.get("id") or s.get("name") or sid)
            pos = (s.get("position") or {})
            x_um = _safe_float(pos.get("x_um")); y_um = _safe_float(pos.get("y_um"))
            spot_ids.append(spot_id)
            spot_paths.append(spath)
            coords_xy.append((x_um if x_um is not None else np.nan,
                              y_um if y_um is not None else np.nan))

        # neighbours (if requested)
        neigh = None
        if spatial_agg_radius is not None or spatial_agg_k is not None:
            neigh = _compute_neighbors_xy(coords_xy, radius=spatial_agg_radius, k=spatial_agg_k)

            # Export adjacency + spot index for potential CAR downstream
            A = _adjacency_from_neigh(neigh, symmetrize=True, self_loops=False)
            if csr_matrix is not None and isinstance(A, csr_matrix):
                save_npz(spatial_dir / f"adjacency_env_{env_id}.npz", A)
            else:
                # Dense fallback: write JSON for portability
                (spatial_dir / f"adjacency_env_{env_id}.json").write_text(json.dumps(A.tolist()))
            (spatial_dir / f"spot_index_env_{env_id}.json").write_text(
                json.dumps({sid: i for i, sid in enumerate(spot_ids)}, indent=2)
            )

        # Precompute neighbourhood metabolite means (env level)
        nn_met_means_by_spot: Dict[str, Dict[str, float]] = {}
        if neigh is not None:
            met_ids: List[str] = []
            met_rows: List[Dict[str, float]] = []
            for spath in spot_paths:
                sobj = _load_spot(spath) or {}
                s = sobj.get("spot") or sobj
                mets = ((s.get("measurements") or {}).get("metabolites") or {}).get("values") or {}
                if not met_ids:
                    met_ids = sorted(map(str, mets.keys()))
                met_rows.append({cid: _safe_float(mets.get(cid)) for cid in met_ids})
            tag = f"r{int(spatial_agg_radius)}" if spatial_agg_radius is not None else f"k{int(spatial_agg_k)}"
            for i, sid in enumerate(spot_ids):
                idxs = neigh[i] if neigh else [i]
                for cid in met_ids:
                    vals = [met_rows[ii].get(cid) for ii in idxs]
                    vals = [v for v in vals if v is not None and np.isfinite(v)]
                    if vals:
                        nn_met_means_by_spot.setdefault(sid, {})[f"nn:met:{cid}_mean_{tag}"] = float(np.mean(vals))

        # Emit rows per (spot, microbe)
        for i, (sid, spath) in enumerate(zip(spot_ids, spot_paths)):
            sobj = _load_spot(spath) or {}
            s = sobj.get("spot") or sobj
            measures = (s.get("measurements") or {})
            microbes = (measures.get("microbes") or {}).get("values") or {}

            # base columns
            pos = (s.get("position") or {})
            base = {
                "spot_id": sid,
                "env_id": env_id,
                "x_um": _safe_float(pos.get("x_um")),
                "y_um": _safe_float(pos.get("y_um")),
                "z_um": _safe_float(pos.get("z_um")),
            }
            # add any top-level annotations from environment if you carry them in your demo (e.g., treatment/cage)
            # (left as-is to avoid changing your structure)

            # coordinate features (optional)
            coord_feats: Dict[str, Any] = {}
            if include_coords:
                for k in ["x_um", "y_um", "z_um"]:
                    if base.get(k) is not None:
                        coord_feats[k] = base[k]
                if coords_poly and (base.get("x_um") is not None) and (base.get("y_um") is not None):
                    x, y = base["x_um"], base["y_um"]
                    coord_feats["x_um2"] = x * x
                    coord_feats["y_um2"] = y * y
                    coord_feats["xy_um2"] = x * y

            for mid, ab in (microbes or {}).items():
                row = dict(base)
                row["microbe"] = str(mid)
                row["abundance"] = _safe_float(ab)

                # coords
                row.update(coord_feats)

                # neighbourhood abundance mean for this microbe
                if neigh is not None:
                    tag = f"r{int(spatial_agg_radius)}" if spatial_agg_radius is not None else f"k{int(spatial_agg_k)}"
                    idxs = neigh[i]
                    nvals = []
                    for ii in idxs:
                        sobj_n = _load_spot(spot_paths[ii]) or {}
                        s_n = sobj_n.get("spot") or sobj_n
                        mic_n = ((s_n.get("measurements") or {}).get("microbes") or {}).get("values") or {}
                        v = _safe_float(mic_n.get(mid))
                        if v is not None and np.isfinite(v):
                            nvals.append(v)
                    if nvals:
                        row[f"nn:abundance_mean_{tag}"] = float(np.mean(nvals))

                # neighbourhood metabolite means
                if nn_met_means_by_spot:
                    for kf, vf in (nn_met_means_by_spot.get(sid, {}) or {}).items():
                        row[kf] = vf

                # target assembly (unchanged capabilities)
                if target:
                    # Look up metabolism for this (spot, microbe)
                    m_spot = (metab_spots.get(sid) or {})
                    m_micro = (m_spot.get("microbes") or {}).get(str(mid)) or {}
                    vals = [_get_target_value(m_micro, tsel) for tsel in target]
                    if len(vals) == 1:
                        row["target"] = vals[0]
                    else:
                        if not target_op:
                            raise typer.BadParameter(
                                "Multiple --target values provided but --target-op not specified "
                                "(valid: sum|mean|min|max|uptake_sum|secretion_sum)."
                            )
                        row["target"] = _combine_targets(vals, target_op)

                rows.append(row)

    # Final table
    df = pd.DataFrame(rows)
    out_table = outdir / "table.csv"
    df.to_csv(out_table, index=False)

    # schema/provenance
    schema = {
        "created_utc": __import__("datetime").datetime.utcnow().isoformat(),
        "system": str(system_yml),
        "metabolism_json": str(metabolism_json) if metabolism_json else None,
        "targets": target or [],
        "target_op": target_op,
        "n_rows": int(df.shape[0]),
        "n_cols": int(df.shape[1]),
        "columns": {c: str(df[c].dtype) for c in df.columns},
        "spatial": {
            "include_coords": bool(include_coords),
            "coords_poly": bool(coords_poly),
            "spatial_agg_radius": spatial_agg_radius,
            "spatial_agg_k": spatial_agg_k,
            "artifacts_dir": str((outdir / "spatial").resolve()),
        },
    }
    (outdir / "schema.json").write_text(json.dumps(schema, indent=2))

    typer.echo(f"Wrote modelling table : {out_table}")
    typer.echo(f"Wrote schema          : {outdir / 'schema.json'}")
    if (spatial_agg_radius is not None or spatial_agg_k is not None):
        typer.echo(f"Spatial artifacts     : {outdir / 'spatial'}")
