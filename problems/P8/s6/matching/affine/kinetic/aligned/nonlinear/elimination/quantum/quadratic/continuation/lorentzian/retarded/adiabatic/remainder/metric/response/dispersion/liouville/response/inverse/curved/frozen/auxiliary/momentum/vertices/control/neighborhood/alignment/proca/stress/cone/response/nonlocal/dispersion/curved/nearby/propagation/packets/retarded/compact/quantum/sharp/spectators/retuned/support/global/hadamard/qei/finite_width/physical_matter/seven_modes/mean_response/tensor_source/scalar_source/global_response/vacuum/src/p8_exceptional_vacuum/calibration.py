"""Exact public calibration domains; no floating point or hidden order choices."""

from fractions import Fraction

import sympy as sp

from . import analytic, family, uniform


def order_point(order=analytic.ORDER):
    if (
        isinstance(order, bool)
        or not isinstance(order, int)
        or order < analytic.ORDER
        or order % 2
    ):
        raise TypeError("Require an even integer switch order n>=1024")
    return {
        "order": order,
        "kappa": sp.Integer(order) / analytic.FIXED_GAMMA,
        "clock_jets_through": order - 1,
        "tensor_domain_lower": -sp.Rational(1, 4 * order),
        "tensor_factor_lower": sp.Rational(1, 2 * order),
        "fixed_gamma": analytic.FIXED_GAMMA,
        "fixed_lambda": analytic.VACUUM_LAMBDA_BAR,
    }


def heavy_remainder_bound(radius):
    if isinstance(radius, bool) or not isinstance(radius, (int, Fraction, sp.Rational)):
        raise TypeError("Require an exact positive rational channel radius")
    radius = sp.Rational(radius)
    D = 2 * analytic.VACUUM_LAMBDA_BAR / analytic.FIXED_GAMMA
    if not 0 < radius < D:
        raise ValueError("The radius must be strictly inside the massive channel gap")
    g2 = 8 * analytic.VACUUM_LAMBDA_BAR * D
    return 3 * g2 * radius**4 / (4 * D**3 * (1 - radius / D))


def bad_cases():
    bad = (True, False, 1.0, sp.Float(1), "1024", sp.I, sp.oo, -sp.oo, sp.zoo, sp.nan)
    calls = [
        ("order_" + str(i), order_point, (v,))
        for i, v in enumerate(bad + (0, -2, 1022, 1025))
    ]
    calls += [
        ("branch_" + str(i), family.constant_branch, (v,))
        for i, v in enumerate(bad + (-1, 2))
    ]
    calls += [
        ("derivative_" + str(i), analytic.derivative_upper, (v,))
        for i, v in enumerate(bad + (-1, 7))
    ]
    D = 2 * analytic.VACUUM_LAMBDA_BAR / analytic.FIXED_GAMMA
    calls += [
        ("radius_" + str(i), heavy_remainder_bound, (v,))
        for i, v in enumerate(bad + (0, -1, D, D + 1))
    ]
    calls += [
        ("rational_envelope_" + str(i), uniform.localized_rational_norm, (v,))
        for i, v in enumerate(
            (
                True,
                1.0,
                sp.Float(1),
                sp.oo,
                sp.nan,
                sp.exp(family.u),
                1 / family.u,
                1 / (1 - family.u**2),
                sp.Symbol("unsupported"),
            )
        )
    ]
    return calls
