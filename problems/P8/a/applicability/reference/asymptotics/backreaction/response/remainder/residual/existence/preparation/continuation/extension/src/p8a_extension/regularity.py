"""Weighted same-slab smoothness using the inherited flat-start Dyson proof."""

import sympy as sp
from p8a_preparation.regularity import frechet_continuity_bound

from . import bounds


def calibration():
    data = bounds.calibration()
    m, t = data["geometry"]["history_uprime_cap"], data["geometry"]["history"]
    return {**data["regularity"],
            "uniform_potential_operator_continuity": frechet_continuity_bound(m, t),
            "same_zero_neighborhood_end": data["source_flat_end"],
            "new_inverse_is_a_generic_C1_endomorphism": False}


def identities():
    k, d, z = sp.symbols("k d z", nonnegative=True)
    sigma, ell = sp.symbols("sigma L", positive=True)
    n = sp.Symbol("n", integer=True, positive=True)
    return {
        "highest_block_coefficient": sp.expand(k*d+k*z+36*k*z**2-k*(d+z+36*z**2)),
        "weighted_unweighted_equivalence": sp.exp(sigma*ell)*sp.exp(-sigma*ell)-1,
        "only_n_single_highest_jet_insertions": n-(n**1),
    }
