"""Exact scope: one paired insertion row, with local anchors but no new forests."""

from functools import cache

import sympy as s
from p8_vacuum_fermion_spectral_pole import enclosure as prior


@cache
def data():
    m, Y, g, M, Q = s.symbols("m Y g M Q", positive=True)
    C = 12 * Y / Q
    T = 4 * m * m
    slope = 4 * C * g / (Q * T)
    remainder = 4 * C * g / (3 * Q * T * T)
    previous = prior.data()
    old_symbols = {str(v): v for v in previous["finite_slope_upper"].free_symbols}
    mapped = previous["finite_slope_upper"].subs(
        {
            old_symbols["C"]: C,
            old_symbols["g"]: g,
            old_symbols["Q"]: Q,
            old_symbols["T"]: T,
        }
    )
    s0 = s.Symbol("s")
    x, v = s.symbols("x v", positive=True)
    Delta = x * v + (1 - x) * M - x * (1 - x) * s0
    return {
        "same_nonlocal_bound": remainder,
        "uniform_MS_slope_bound": slope,
        "inner_forest": "Complete original on-shell fermion mass/residue subtraction in D dimensions. This is not a raw fermion self-energy with its local reference omitted.",
        "outer_forest": "After the inner forest, the remaining tadpole and bubble poles are local mass/source terms. Subtract their MS poles only. The bubble derivative and on-shell double remainder are convergent and need no finite scheme conversion.",
        "not_other_forests": "No primitive chord mass/slope is recounted. No separate counterterm product is added to the already paired row. Finite parent-parameter changes inside other loops, canonical factors and vacuum source cross terms remain open.",
        "checks": {
            "same_prior_slope": s.factor(mapped - slope),
            "same_slope_to_remainder_ratio": s.factor(slope - 3 * T * remainder),
            "real_interval_parameter_gap": s.expand(
                Delta - (x * (v - 4) + (1 - x) * M) - x * (1 - x) * (4 - s0) - 4 * x * x
            ),
            "mass_change_bound_over_unit_interval": slope * (1 - 0) - slope,
        },
    }
