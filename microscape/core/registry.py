
from typing import Callable, Dict
KineticsStep = Callable[[dict, dict, float, float], dict]
_REGISTRY: Dict[str, KineticsStep] = {}
def register(name: str):
    def deco(fn: KineticsStep):
        _REGISTRY[name] = fn
        return fn
    return deco
def get(name: str) -> KineticsStep:
    if name not in _REGISTRY:
        raise KeyError(f"Unknown engine '{name}'. Available: {', '.join(sorted(_REGISTRY))}")
    return _REGISTRY[name]
