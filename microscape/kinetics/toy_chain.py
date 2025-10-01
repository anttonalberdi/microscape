
from ..core.registry import register
@register("toy")
def step(fields: dict, cfg: dict, dt: float, alpha: float) -> dict:
    glc = fields.get("glc"); lac = fields.get("lac"); ac = fields.get("ac"); but = fields.get("but")
    if glc is None or lac is None or ac is None or but is None:
        raise ValueError("toy engine expects fields: glc, lac, ac, but")
    uptake_k = (cfg.get("constraints") or {}).get("uptake_k") or {}
    k_glc = float(uptake_k.get("glc", 0.05))
    k_lac = float(uptake_k.get("lac", 0.04))
    k_ac  = float(uptake_k.get("ac",  0.03))
    y_fd_lac, y_fd_ac = 0.6, 0.4
    y_lu_but, y_bp_but = 0.9, 0.8
    r_fd = alpha * k_glc * glc
    r_lu = alpha * k_lac * lac
    r_bp = alpha * k_ac  * ac
    return {"glc": -r_fd, "lac": y_fd_lac*r_fd - r_lu, "ac": y_fd_ac*r_fd - r_bp, "but": y_lu_but*r_lu + y_bp_but*r_bp}
