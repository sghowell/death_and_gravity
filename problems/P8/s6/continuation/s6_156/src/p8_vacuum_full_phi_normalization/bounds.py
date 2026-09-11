"""Exact upper estimates for the complete second normalization and pole."""

from functools import cache

import sympy as s
from p8_vacuum_two_loop_finite_contact.calibration import rational


def normalization_upper(hybrid_slope, k0, r0):
    vals = tuple(map(rational, (hybrid_slope, k0, r0)))
    if min(vals) < 0:
        raise ValueError("Need nonnegative exact slope and pole-coefficient bounds")
    T, k, R = vals
    return T + k * R


@cache
def data():
    T, k, R, Bs, Bf, B2 = s.symbols("T k R Bscalar Bferm Btwo", nonnegative=True)
    full = T + k * R
    pole = Bs + Bf + B2 + k * (2 * Bs + Bf)
    return {
        "complete_MS_second_normalization_upper": full,
        "finite_parameter_correction_to_OS_coefficient_upper": k * (2 * Bs + Bf),
        "full_one_plus_two_loop_OS_coefficient_upper": pole,
        "checks": {
            "all_second_normalization_terms_retained": full - T - k * R,
            "both_first_loop_parameter_directions": s.expand(
                k * (2 * Bs + Bf) - 2 * k * Bs - k * Bf
            ),
            "canonical_pole_sum_with_conversion": pole
            - Bs
            - Bf
            - B2
            - k * (2 * Bs + Bf),
            "normalization_lower_gap": (1 - k - full) + k + full - 1,
        },
        "scope": "Only nonnegative absolute enclosures are combined, so differing inherited self-energy sign conventions cannot cancel an error. The finite MS-to-hybrid re-expansion and every complete quadratic family are included. No physical higher-loop remainder is inferred from a finite-order polynomial normalization.",
    }
