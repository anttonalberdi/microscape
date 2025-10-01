
from __future__ import annotations
from typing import Dict, Tuple, Callable
import numpy as np
from ..io.graph_config import parse_graph, load_graph_yaml
from microscape.space import build_laplacians, explicit_diffusion_step, stable_dt_upper_bound

KineticsStep = Callable[[Dict[str, np.ndarray], dict, float, float], Dict[str, np.ndarray]]

def simulate_graph(config_path: str, kinetics_step: KineticsStep, *, progress=None) -> Tuple[Dict[str, np.ndarray], dict]:
    cfg = load_graph_yaml(config_path)
    g = parse_graph(cfg)
    field_names = sorted({k for nd in g.nodes for k in nd.init.keys()})
    N = len(g.nodes)
    fields = {f: np.zeros(N, dtype=float) for f in field_names}
    for idx, nd in enumerate(g.nodes):
        for k,v in nd.init.items():
            fields[k][idx] = float(v)
    Ls = build_laplacians(N, g.edges, g.diffusion_um2_s, field_names)
    steps = int(cfg.get("runtime", {}).get("steps", 200))
    dt = float(cfg.get("runtime", {}).get("dt_s", 2.0))
    alpha = float(cfg.get("constraints", {}).get("alpha_flux_to_dC", 5e-3))
    for f in field_names:
        Dt = float(g.diffusion_um2_s.get(f, 0.0))
        dt = min(dt, stable_dt_upper_bound(Ls[f], Dt, safety=0.45))
    for t in range(steps):
        dC = kinetics_step(fields, cfg, dt, alpha)
        for f in field_names:
            decay = float(g.decay_s.get(f, 0.0))
            fields[f] = np.maximum(0.0, fields[f] + dt * (dC[f] - decay * fields[f]))
        explicit_diffusion_step(fields, Ls, g.diffusion_um2_s, dt)
        if progress and (t % max(1, steps//100) == 0):
            progress(int(100*(t+1)/steps), 100)
    summary = {f"mean_{k}": float(v.mean()) for k, v in fields.items()}
    summary.update({"n_steps": float(steps), "n_nodes": float(N)})
    return fields, summary
