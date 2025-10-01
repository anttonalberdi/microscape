
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import yaml

@dataclass
class Node:
    id: str
    pos: Optional[Tuple[float, float, float]]
    volume_nl: float
    init: Dict[str, float]
    guilds: Optional[Dict[str, float]] = None
    transcripts: Optional[Dict[str, float]] = None
    region: Optional[str] = None

@dataclass
class Edge:
    i: int
    j: int
    weight: float
    D_scale: Optional[Dict[str, float]] = None
    bias: Optional[Dict[str, float]] = None

@dataclass
class GraphSpace:
    nodes: List[Node]
    edges: List[Edge]
    diffusion_um2_s: Dict[str, float]
    decay_s: Dict[str, float]

def load_graph_yaml(path: str | Path) -> dict:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Graph config not found: {path}")
    return yaml.safe_load(p.read_text())

def parse_graph(cfg: dict) -> GraphSpace:
    space = cfg.get("space", {})
    if space.get("type") != "graph":
        raise ValueError("space.type must be 'graph'")
    # nodes
    nodes: List[Node] = []
    id_to_idx: Dict[str, int] = {}
    for idx, nd in enumerate(space.get("nodes", [])):
        id_to_idx[nd["id"]] = idx
        pos = nd.get("pos_um")
        if pos is not None:
            pos = tuple(float(x) for x in pos)
            if len(pos) == 2:
                pos = (pos[0], pos[1], 0.0)
        nodes.append(Node(
            id=nd["id"],
            pos=pos,
            volume_nl=float(nd.get("volume_nl", 0.05)),
            init={k: float(v) for k,v in (nd.get("init") or {}).items()},
            guilds=nd.get("guilds"),
            transcripts=nd.get("transcripts"),
            region=nd.get("region"),
        ))
    # edges
    edges: List[Edge] = []
    for e in space.get("edges", []):
        i = id_to_idx[e["i"]]; j = id_to_idx[e["j"]]
        edges.append(Edge(
            i=i, j=j,
            weight=float(e.get("weight", 1.0)),
            D_scale=e.get("D_scale"),
            bias=e.get("bias"),
        ))
    diffusion = space.get("diffusion_um2_s", {})
    decay = space.get("decay_s", {})
    return GraphSpace(nodes, edges, diffusion, decay)
