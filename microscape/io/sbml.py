from __future__ import annotations
from pathlib import Path
from typing import Dict, Any

def require_cobra():
    try:
        import cobra  # noqa: F401
    except Exception as e:
        raise RuntimeError(
            "COBRApy is required for SBML demos. Install with:\n"
            "  pip install cobra optlang\n"
            "or via conda:\n"
            "  conda install -c conda-forge cobra"
        ) from e

def load_sbml_model(path: str):
    """Load an SBML model with COBRApy, returning the model object."""
    require_cobra()
    from cobra.io import read_sbml_model
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"SBML model not found: {p}")
    return read_sbml_model(str(p))

def map_exchanges(model, wanted: Dict[str, str]) -> Dict[str, Any]:
    """Return a dict mapping field_name -> cobra Reaction for exchanges.
    wanted: e.g., {'glc': 'EX_glc__D_e', 'lac': 'EX_lac__L_e'}
    """
    rxns = {}
    for field, rxn_id in wanted.items():
        try:
            rxns[field] = model.reactions.get_by_id(rxn_id)
        except KeyError:
            raise KeyError(f"Exchange '{rxn_id}' not found in model '{model.id}' for field '{field}'.")
    return rxns
