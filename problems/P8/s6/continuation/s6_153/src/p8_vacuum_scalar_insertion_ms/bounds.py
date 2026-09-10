"""Positive integral control at the actual fixed MS scale mu=mF."""

from functools import cache

import sympy as s
from p8_vacuum_two_loop_finite_contact.calibration import rational


def enclosure(L, g, M, tree_b2, Q=144):
    L, g, M, tree_b2, Q = map(rational, (L, g, M, tree_b2, Q))
    if (
        L <= 0
        or g <= 0
        or tree_b2 <= 0
        or not 32 <= M <= 10**400
        or not 0 < Q <= 144
        or not g / M < L / 3
    ):
        raise ValueError("Need L,g,tree>0; 32<=M<=10^400; g/M<L/3; 0<Q_lower<=144")
    F = g * 1200 / (Q**2 * M)
    return {
        "alpha_zero_upper": g / (2 * Q * M),
        "finite_F_alpha_upper": F,
        "local_conversion_b2_relative_upper": 2 * L * F,
        "local_conversion_b2_absolute_upper": tree_b2 * 2 * L * F,
    }


@cache
def data():
    x, M, g, Q, ell = s.symbols("x M g Q ell", positive=True)
    D = x * M + (1 - x) ** 2
    b = x * (1 - x) / D
    pref = g / Q**2
    return {
        "finite_reference_integrand": pref * b * (2 * ell - s.log(D)),
        "finite_reference_upper": g * ell / (Q**2 * M),
        "positive_scale_hierarchy": "mu^2=mF^2=10^400 >= M >= 32",
        "checks": {
            "denominator_monotone_derivative": s.diff(D, x) - (M - 2 + 2 * x),
            "finite_weight_upper_integral": pref
            * 2
            * ell
            * s.integrate((1 - x) / M, (x, 0, 1))
            - g * ell / (Q**2 * M),
            "forward_ratio_follows_from_actual_domain_margin": s.expand(
                3 * (M - 2) - 2 * M - (M - 6)
            ),
            "two_loop_measure_power": 144**2 - 20736,
            "scale_logarithm_bound_prefactor": 2 * 200 * 3 - 1200,
        },
        "scope": "For the fixed scale, 1<=D<=M<=mu^2 implies ell<=2ell-logD<=2ell, so the finite reference is positive. The bound includes the slope's epsilon coefficient. It is a fixed-order local-reference estimate, not a physical higher-loop remainder.",
    }
