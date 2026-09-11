"""Private exact mixed metric-jet seminorm bounds for the fixed shear family."""

from functools import cache
from math import comb, factorial

import sympy as s
from p8_vacuum_affine_matrix_adiabatic import jets as prior

TIME = 4
PARAMETERS = 3
G = s.Integer(2) * 10**5


def exponential(j, a):
    return 2 * sum(
        comb(j, l) * s.Integer(a) ** (j - l) * prior.bell(l) for l in range(j + 1)
    )


def scale(j, inverse=False):
    if inverse:
        return prior.scale_bound(2, j)
    t = s.Symbol("t", real=True)
    polynomial = s.Poly(s.diff((1 + t * t) ** 4, t, j), t)
    return sum(
        abs(coefficient) * s.Rational(1, 2) ** power[0]
        for power, coefficient in polynomial.terms()
    )


@cache
def constants():
    raw = {}
    shifted = {}
    for inverse in (False, True):
        for j in range(TIME + 1):
            for a in range(PARAMETERS + 1):
                raw[inverse, j, a] = sum(
                    comb(j, l) * scale(j - l, inverse) * exponential(l, a)
                    for l in range(j + 1)
                )
        for r in range(TIME + 1):
            value = sum(
                6 * raw[inverse, j + r, a] / (factorial(j) * factorial(a))
                for j in range(TIME - r + 1)
                for a in range(PARAMETERS + 1)
            )
            if r == 0:
                value += 1
            shifted[inverse, r] = s.factor(value)
    detector = 30 * s.Rational(25, 16) ** 2
    return {
        "absolute_mixed_exponential": {
            (j, a): exponential(j, a)
            for j in range(TIME + 1)
            for a in range(PARAMETERS + 1)
        },
        "raw_metric_operator_bounds": raw,
        "all_time_shifted_component_jet_bounds": shifted,
        "detector_metric_variation_C0_2_component_bound": detector,
    }


@cache
def data():
    c = constants()
    t = s.Symbol("t", real=True)
    checks = {}
    for j in range(5):
        polynomial = s.diff((1 + t * t) ** 4, t, j)
        checks[f"positive_polynomial_metric_scale_derivative_{j}"] = s.factor(
            polynomial.subs(t, s.Rational(1, 2)) - scale(j)
        )
    checks["absolute_zero_jet_exponential_is_not_relative_one"] = exponential(0, 0) - 2
    return {
        "metric_jet_display": G,
        "constants": c,
        "checks": checks,
        "gates": {
            "all_metric_and_inverse_time_shifted_jets": all(
                v < G for v in c["all_time_shifted_component_jet_bounds"].values()
            ),
            "detector_no_time_derivative_required": c[
                "detector_metric_variation_C0_2_component_bound"
            ]
            < G,
            "conservative_tensor_algebra_domain": G >= 3,
        },
    }
