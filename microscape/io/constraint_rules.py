from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional, Any
import yaml

@dataclass
class NormalizeSpec:
    mode: str = "per_microbe_max"     # per_microbe_max | per_model_max | constant
    constant: float = 1000.0          # used if mode == constant

@dataclass
class GPRActivitySpec:
    logic: str = "and_min_or_max"     # only option we implement now
    normalize: NormalizeSpec = field(default_factory=NormalizeSpec)
    threshold_TPM: float = 1.0
    floor_activity: float = 0.0
    cap_activity: float = 1.0

@dataclass
class BoundScalingSpec:
    min_irrev_ub: float = 0.0
    min_rev_mag: float = 0.0

@dataclass
class UptakeSpec:
    enabled: bool = True
    mM_per_uptake: float = 1.0
    max_uptake: float = 10.0
    secretion_upper: float = 1000.0
    metabolite_map: Dict[str, str] = field(default_factory=dict)

@dataclass
class WriteModelsSpec:
    enabled: bool = True
    dir: str = "constrained_models"

@dataclass
class ConstraintRules:
    solver: str = "glpk"
    gpr_activity: GPRActivitySpec = field(default_factory=GPRActivitySpec)
    bound_scaling: BoundScalingSpec = field(default_factory=BoundScalingSpec)
    uptake: UptakeSpec = field(default_factory=UptakeSpec)
    write_models: WriteModelsSpec = field(default_factory=WriteModelsSpec)

def load_constraint_rules(yml: Path) -> ConstraintRules:
    data = yaml.safe_load(yml.read_text())
    d = (data or {}).get("constraints", {})
    # parse nested
    def _norm(dn: Optional[dict]) -> NormalizeSpec:
        dn = dn or {}
        return NormalizeSpec(mode=str(dn.get("mode","per_microbe_max")),
                             constant=float(dn.get("constant", 1000.0)))
    def _gpr(dg: Optional[dict]) -> GPRActivitySpec:
        dg = dg or {}
        ns = _norm(dg.get("normalize"))
        return GPRActivitySpec(
            logic=str(dg.get("logic","and_min_or_max")),
            normalize=ns,
            threshold_TPM=float(dg.get("threshold_TPM",1.0)),
            floor_activity=float(dg.get("floor_activity",0.0)),
            cap_activity=float(dg.get("cap_activity",1.0)),
        )
    def _bs(db: Optional[dict]) -> BoundScalingSpec:
        db = db or {}
        return BoundScalingSpec(
            min_irrev_ub=float(db.get("min_irrev_ub", 0.0)),
            min_rev_mag=float(db.get("min_rev_mag", 0.0))
        )
    def _up(du: Optional[dict]) -> UptakeSpec:
        du = du or {}
        return UptakeSpec(
            enabled=bool(du.get("enabled", True)),
            mM_per_uptake=float(du.get("mM_per_uptake", 1.0)),
            max_uptake=float(du.get("max_uptake", 10.0)),
            secretion_upper=float(du.get("secretion_upper", 1000.0)),
            metabolite_map=dict(du.get("metabolite_map", {})),
        )
    def _wm(dw: Optional[dict]) -> WriteModelsSpec:
        dw = dw or {}
        return WriteModelsSpec(
            enabled=bool(dw.get("enabled", True)),
            dir=str(dw.get("dir", "constrained_models"))
        )
    return ConstraintRules(
        solver=str(d.get("solver","glpk")),
        gpr_activity=_gpr(d.get("gpr_activity")),
        bound_scaling=_bs(d.get("bound_scaling")),
        uptake=_up(d.get("uptake")),
        write_models=_wm(d.get("write_models")),
    )
