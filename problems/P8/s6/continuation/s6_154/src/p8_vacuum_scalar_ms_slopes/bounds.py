"""Exact rational upper enclosures at the fixed mu^2=10^400 reference."""

from functools import cache

import sympy as s
from p8_vacuum_two_loop_finite_contact.calibration import rational


def enclosure(L, g, M, Q=144):
    L, g, M, Q = map(rational, (L, g, M, Q))
    if L <= 0 or g <= 0 or not 32 < M <= 10**400 or not 0 < Q <= 144:
        raise ValueError("Need L,g>0; 32<M<=10^400; 0<Q_lower<=144")
    ell = s.Integer(1200)
    alpha = g / (2 * Q * M)
    groups = {
        "local_sunset_MS": L**2 / (6 * Q**2) * (ell / 2 + 3 + s.Rational(1, 16)),
        "mixed_sunset_anchored": s.Rational(33, 2) * L * g / Q**2,
        "same_heavy_sunset": s.Rational(5, 4) * g**2 / Q**2,
        "different_heavy_sunset": 28 * g**2 / Q**2,
        "nested_decaying_remainder": 2 * g**2 / Q**2,
        "nested_constant_alpha": alpha**2,
        "proper_interaction_MS_conversion": L * g * ell / (2 * Q**2 * M)
        + g**2 * ell / (4 * Q**2 * M**2),
    }
    return {
        "finite_MS_slope_group_uppers": groups,
        "finite_MS_slope_absolute_upper": sum(groups.values()),
        "additional_OS_quadratic_upper": ell
        * (L * g + g**2 / M)
        / (6 * Q**2 * M**2 * (1 - 1 / M)),
        "same_inner_alpha_upper": alpha,
    }


@cache
def data():
    r, t = s.symbols("r t", positive=True)
    v = 1 + t + r * t
    W = t / v**3
    return {
        "local_T1_absolute_upper": s.Integer(6),
        "local_C_positive_upper": s.Rational(1, 16),
        "scale_logarithm_upper": s.Integer(1200),
        "checks": {
            "log_r_integrable_majorant": 6
            * s.integrate(t, (t, 0, 1))
            * s.integrate(-s.log(r), (r, 0, 1))
            - 3,
            "bounded_sector_logs_allowance": 3 * 2 * s.Rational(1, 2) - 3,
            "log_one_minus_P_upper": s.Rational(1, 2)
            * s.Rational(1, 9)
            / (1 - s.Rational(1, 9))
            - s.Rational(1, 16),
            "local_weight_below_t_gap": s.factor(t - W - t * (v**3 - 1) / v**3),
            "scale_log_bound": 400 * 3 - 1200,
            "two_loop_measure_lower": 144**2 - 20736,
        },
        "scope": "Positive simplex/sector bounds: U>0, P<=1/9, integral w=1/2, |T1|<6 and 0<C<=1/16. The large heavy mass is dropped only in already-convergent inherited estimates. These unequal conservative bounds do not diagnose perturbative growth.",
    }
