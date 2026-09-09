"""Exact complete-geometry negative control: caps alone do not force an endpoint."""

from functools import cache

import sympy as sp

from . import costs


@cache
def data():
    x = sp.Symbol("x", real=True)
    K = sp.Integer(100)
    q = (
        2 * x
        + (2 * K + 4) * x * x
        + (K * K + 4 * K + sp.Rational(40, 3)) * x**3
        + (K**3 / 3 + 2 * K * K + sp.Rational(40, 3) * K + sp.Rational(160, 3)) * x**4
    )
    A = 1 - sp.exp(-K * x) * q
    radiation = (1 - 6 * x) ** sp.Rational(1, 3)
    checks = {
        f"complete_control_matches_radiation_affine_jet_{j}": sp.diff(
            A - radiation, x, j
        ).subs(x, 0)
        for j in range(5)
    }
    y = sp.Symbol("y", nonnegative=True)
    scaled = sp.Poly(K * q.subs(x, y / K), y)
    S = sum(coeff * sp.factorial(power[0]) for power, coeff in scaled.terms())
    future = {j: sp.factor(2**j * K ** (j - 1) * S) for j in range(1, 5)}
    past = {}
    for j in range(1, 5):
        p = sp.Poly(sp.expand(sp.diff(A, x, j) * sp.exp(K * x)), x)
        past[j] = sp.Rational(10, 9) * sum(
            abs(coeff) * costs.R ** power[0] for power, coeff in p.terms()
        )
    past_value = sp.Rational(10, 9) * sum(
        abs(coeff) * costs.R ** power[0] for power, coeff in sp.Poly(q, x).terms()
    )
    return {
        "A": A,
        "future_S": S,
        "future_A_lower": 1 - S / K,
        "future_derivative_bounds": future,
        "past_absolute_A_minus_one_bound": past_value,
        "past_derivative_bounds": past,
        "past_derivative_margins": {
            j: costs.PAST[2][j - 1] - past[j] for j in range(1, 5)
        },
        "future_derivative_margins": {
            j: costs.FUTURE[2][j - 1] - future[j] for j in range(1, 5)
        },
        "checks": checks,
    }
