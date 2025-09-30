from __future__ import annotations
from typing import Dict

def run_fba_with_bounds(model, exchange_bounds: Dict[str, float], use_pfba: bool = True):
    """Run (p)FBA after setting exchange lower bounds from exchange_bounds.
    exchange_bounds: mapping reaction_id -> uptake cap (positive number); will set lower_bound = -cap
    Returns a cobra Solution.
    """
    from contextlib import contextmanager
    from cobra.flux_analysis import pfba as cobra_pfba

    @contextmanager
    def temp_bounds(m, bounds: Dict[str, float]):
        with m:
            for rxn_id, cap in bounds.items():
                try:
                    rxn = m.reactions.get_by_id(rxn_id)
                except KeyError:
                    continue
                # BiGG convention: uptake is negative flux
                rxn.lower_bound = -abs(cap)
            yield m

    with temp_bounds(model, exchange_bounds):
        if use_pfba:
            sol = cobra_pfba(model)
        else:
            sol = model.optimize()
    return sol
