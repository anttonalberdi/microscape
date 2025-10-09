from __future__ import annotations
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
import json
import math
import numpy as np
import pandas as pd
import typer
import yaml
from scipy.sparse import csr_matrix, save_npz  # falls back to local if SciPy missing (see below)

app = typer.Typer(add_completion=False, no_args_is_help=True)

# ---- Helpers ----------------------------------------------------------------

def _safe_float(x):
    try:
        return float(x)
    except Exception:
        return None

def _iter_env_files(system_yml: Path) -> List[Path]:
    doc = yaml.safe_load(Path(system_yml).read_text())
    sysn = (doc.get("system") or {})
    envs = sysn.get("environments") or []
    base = Path(system_yml).parent
    return [ (base / e).resolve() for e in envs ]

def _iter_spots_for_env(env_file: Path, system_root: Path) -> List[Tuple[str, Path]]:
    env = yaml.safe_load(Path(env_file).read_text()) or {}
    e = env.get("environment") or {}
    spots = e.get("spots") or []
    base = system_root
    out = []
    for s in spots:
        # allow both direct path or id->path
        if isinstance(s, str):
            sid = Path(s).stem
            sp = (base / "spots" / s).resolve() if not s.endswith(".yml") else (base / s).resolve()
        else:
            sid = s.get("id") or s.get("name")
            sp = (base / s.get("path")).resolve()
        out.append((str(sid), sp))
    return out

def _load_spot(spot_file: Path) -> Dict[str, Any]:
    return yaml.safe_load(Path(spot_file).read_text()) or {}

def _compute_neighbors_xy(coords: List[Tuple[float,float]], radius: Optional[float]=None, k: Optional[int]=None) -> List[List[int]]:
    """Neighbour lists for each point. If radius, use ball; if k, use kNN; if neither, self only."""
    n = len(coords)
    if n == 0:
        return []
    xs = np.array([c[0] for c in coords], dtype=float)
    ys = np.array([c[1] for c in coords], dtype=float)
    valid = np.isfinite(xs) & np.isfinite(ys)
    # brute force to keep dependencies minimal; n (spots per env) is typically small/medium
    neigh: List[List[int]] = [[] for _ in range(n)]
    if radius is None and (k is None or k <= 1):
        for i in range(n):
            neigh[i] = [i]
        return neigh
    # fill distances
    for i in range(n):
        if not valid[i]:
            neigh[i] = []
            continue
        di = []
        xi, yi = xs[i], ys[i]
        for j in range(n):
            if not valid[j]:
                continue
            dx, dy = xi - xs[j], yi - ys[j]
            d = math.hypot(dx, dy)
            if radius is not None:
                if d <= float(radius):
                    di.append(j)
            else:
                di.append((d, j))
        if k is not None:
            di = [j for _, j in sorted(di)[:max(1, int(k))]]
        neigh[i] = di if radius is not None else di
    return neigh

def _adjacency_from_neigh(neigh: List[List[int]], symmetrize: bool=True, self_loops: bool=False) -> csr_matrix:
    n = len(neigh)
    rows, cols, data = [], [], []
    for i, js in enumerate(neigh):
        for j in js:
            if (not self_loops) and (i == j):
                continue
            rows.append(i); cols.append(j); data.append(1)
    A = csr_matrix((data, (rows, cols)), shape=(n, n))
    if symmetrize:
        A = ((A + A.T) > 0).astype(int)
    return A

# ---- CLI --------------------------------------------------------------------

@app.command("build")
def build_cmd(
    system_yml: Path = typer.Argument(..., help="Path to system.yml"),
    outdir: Path = typer.Option("outputs/model", help="Directory for outputs."),
    metabolism_json: Optional[Path] = typer.Option(None, "--metabolism-json", help="(Optional) metabolism JSON for selecting targets elsewhere."),
    target: Optional[str] = typer.Option(None, help="(Handled elsewhere in your pipeline)"),
    include_coords: bool = typer.Option(False, help="Include x_um,y_um,z_um."),
    coords_poly: bool = typer.Option(False, help="Add x_um^2,y_um^2, x*y."),
    spatial_agg_radius: Optional[float] = typer.Option(None, help="Radius (µm) for neighbourhood mean features."),
    spatial_agg_k: Optional[int] = typer.Option(None, help="k for kNN neighbourhood mean features."),
):
    """
    Assemble the modeling table. This version can optionally add:
      - raw coordinates (x_um,y_um,z_um) and simple polynomials,
      - neighbourhood mean features using radius or kNN per environment,
      - per-environment adjacency & spot index for spatial CAR models.
    """
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    spatial_dir = outdir / "spatial"
    spatial_dir.mkdir(exist_ok=True)

    system_root = Path(system_yml).parent.resolve()
    env_files = _iter_env_files(system_yml)
    rows: List[Dict[str, Any]] = []

    # Iterate environments
    for env_file in env_files:
        env_doc = yaml.safe_load(Path(env_file).read_text()) or {}
        env_meta = (env_doc.get("environment") or {})
        env_id = str(env_meta.get("id") or Path(env_file).stem)

        # Gather all spots in this env (order matters for adjacency)
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

        # Build neighbour structure (if requested)
        neigh = None
        if spatial_agg_radius is not None or spatial_agg_k is not None:
            neigh = _compute_neighbors_xy(coords_xy, radius=spatial_agg_radius, k=spatial_agg_k)

        # Export adjacency + spot index for CAR
        if neigh is not None:
            A = _adjacency_from_neigh(neigh, symmetrize=True, self_loops=False)
            save_npz(spatial_dir / f"adjacency_env_{env_id}.npz", A)
            (spatial_dir / f"spot_index_env_{env_id}.json").write_text(
                json.dumps({sid: i for i, sid in enumerate(spot_ids)}, indent=2)
            )

        # Precompute neighbourhood means for metabolites (optional)
        nn_met_means_by_spot: Dict[str, Dict[str, float]] = {}
        if neigh is not None:
            # discover metabolite ID order from the first spot that has them
            met_ids: List[str] = []
            # Build per-spot metabolite dicts
            met_rows: List[Dict[str, float]] = []
            for spath in spot_paths:
                sobj = _load_spot(spath) or {}
                s = sobj.get("spot") or sobj
                mets = ((s.get("measurements") or {}).get("metabolites") or {}).get("values") or {}
                if not met_ids:
                    met_ids = sorted(map(str, mets.keys()))
                met_rows.append({cid: _safe_float(mets.get(cid)) for cid in met_ids})

            # compute NN mean per spot per metabolite
            tag = f"r{int(spatial_agg_radius)}" if spatial_agg_radius is not None else f"k{int(spatial_agg_k)}"
            for i, sid in enumerate(spot_ids):
                idxs = neigh[i] if neigh else [i]
                vals_cache = {}
                for cid in met_ids:
                    vals = [met_rows[ii].get(cid) for ii in idxs]
                    vals = [v for v in vals if v is not None and np.isfinite(v)]
                    if vals:
                        nn_met_means_by_spot.setdefault(sid, {})[f"nn:met:{cid}_mean_{tag}"] = float(np.mean(vals))

        # Now write rows for each (spot, microbe)
        for i, (sid, spath) in enumerate(zip(spot_ids, spot_paths)):
            sobj = _load_spot(spath) or {}
            s = sobj.get("spot") or sobj
            measures = (s.get("measurements") or {})
            microbes = (measures.get("microbes") or {}).get("values") or {}

            # Base columns common to every microbe row
            pos = (s.get("position") or {})
            base = {
                "spot_id": sid,
                "env_id": env_id,
                "x_um": _safe_float(pos.get("x_um")),
                "y_um": _safe_float(pos.get("y_um")),
                "z_um": _safe_float(pos.get("z_um")),
                # You likely already add: treatment, cage, etc. if present in your system
            }

            # Coordinates as features
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

            # For each microbe with an abundance value, emit a row
            for mid, ab in (microbes or {}).items():
                row = dict(base)
                row["microbe"] = str(mid)
                row["abundance"] = _safe_float(ab)

                # add coordinate features
                row.update(coord_feats)

                # neighbourhood mean of ABUNDANCE for the same microbe
                if neigh is not None:
                    tag = f"r{int(spatial_agg_radius)}" if spatial_agg_radius is not None else f"k{int(spatial_agg_k)}"
                    idxs = neigh[i]
                    nvals = []
                    for ii in idxs:
                        sobj_n = _load_spot(spot_paths[ii]) or {}
                        s_n = sobj_n.get("spot") or sobj_n
                        mic_n = ((s_n.get("measurements") or {}).get("microbes") or {}).get("values") or {}
                        v = mic_n.get(mid)
                        v = _safe_float(v)
                        if v is not None and np.isfinite(v):
                            nvals.append(v)
                    if nvals:
                        row[f"nn:abundance_mean_{tag}"] = float(np.mean(nvals))

                # neighbourhood metabolite means (precomputed per spot)
                if nn_met_means_by_spot:
                    for k, v in (nn_met_means_by_spot.get(sid, {}) or {}).items():
                        row[k] = v

                rows.append(row)

    # Final table
    df = pd.DataFrame(rows)
    out_table = outdir / "table.csv"
    df.to_csv(out_table, index=False)

    # Record spatial manifest for later
    manifest = {
        "include_coords": bool(include_coords),
        "coords_poly": bool(coords_poly),
        "spatial_agg_radius": spatial_agg_radius,
        "spatial_agg_k": spatial_agg_k,
        "car_adjacency_dir": str(spatial_dir),
        "environments": [Path(e).stem for e in env_files],
    }
    (outdir / "spatial_manifest.json").write_text(json.dumps(manifest, indent=2))

    typer.echo(f"Wrote modeling table: {out_table}")
    typer.echo(f"Spatial artifacts: {spatial_dir}")
