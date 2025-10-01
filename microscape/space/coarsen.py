from __future__ import annotations
from typing import Dict, List, Tuple
from .graph import Edge

def coarsen(nodes, edges: List[Edge], groups: List[List[int]]):
    """
    Merge voxel indices in 'groups' into super-voxels.
    - Volume: sum
    - Guilds: keep as-is from representative (or implement weighted merge as needed)
    - Edges: sum conductances between merged groups
    """
    # old->new mapping
    mapping: Dict[int, int] = {}
    for new_i, grp in enumerate(groups):
        for old_i in grp:
            mapping[old_i] = new_i

    # Coarsened nodes (minimal: keep id/pos/volume; extend as needed)
    new_nodes = []
    for new_i, grp in enumerate(groups):
        rep = nodes[grp[0]]
        vol = sum(getattr(nodes[i], "volume_nl", 0.0) for i in grp)
        # Construct your project’s Node class; here we reuse the same type
        new_nodes.append(type(rep)(
            id=f"g{new_i}",
            pos=getattr(rep, "pos", None),
            volume_nl=vol,
            init=getattr(rep, "init", {}),
            guilds=getattr(rep, "guilds", {}),
            transcripts=getattr(rep, "transcripts", {}),
            region=getattr(rep, "region", None),
        ))

    # Aggregate edges between groups
    acc: Dict[Tuple[int, int], Dict] = {}
    for e in edges:
        i2, j2 = mapping[e.i], mapping[e.j]
        if i2 == j2:
            continue
        key = (min(i2, j2), max(i2, j2))
        rec = acc.setdefault(key, {"i": key[0], "j": key[1], "weight": 0.0, "D_scale": {}})
        rec["weight"] += float(e.weight)
        for k, v in (e.D_scale or {}).items():
            rec["D_scale"][k] = rec["D_scale"].get(k, 0.0) + float(v)

    new_edges = [Edge(**rec) for rec in acc.values()]
    return new_nodes, new_edges
