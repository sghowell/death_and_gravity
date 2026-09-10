"""Exact order, connected-domain and unrestricted-Hessian calibrations."""

from fractions import Fraction

import sympy as sp
from p8_exceptional_vacuum import calibration as previous

from . import covariance, family


def point(order=family.N, physical_X=1):
    row = previous.order_point(order)
    if isinstance(physical_X, bool) or not isinstance(
        physical_X, (int, Fraction, sp.Rational)
    ):
        raise TypeError("An exact rational kinetic invariant is required")
    X = sp.Rational(physical_X)
    if not -sp.Rational(1, 4 * order) < X < sp.Rational(6, 5):
        raise ValueError("The stated open connected field domain is required")
    return {
        **row,
        "physical_X": X,
        "tensor_factor_lower": sp.Rational(1, 2) + sp.Rational(1, 8 * order),
        "additional_quotient_factor_lower": sp.Rational(1, 4 * order),
        "regular_q_bound_coefficient": 32 * order,
    }


def bad_cases():
    invalid = (True, False, 1.0, sp.Float(1), "1", sp.I, sp.oo, -sp.oo, sp.zoo, sp.nan)
    calls = [
        ("order_" + str(i), point, (v,))
        for i, v in enumerate(invalid + (0, 1022, 1025))
    ]
    calls += [
        ("X_" + str(i), point, (family.N, v))
        for i, v in enumerate(
            invalid + (-1, -sp.Rational(1, 4 * family.N), sp.Rational(6, 5), 2)
        )
    ]
    calls += [
        ("gradient_shape_" + str(i), covariance.hessian, (v, 1, 0, 0))
        for i, v in enumerate((None, "0000", (), (0, 0, 0), (0, 0, 0, 0, 0)))
    ]
    for i, v in enumerate(invalid):
        calls.append(
            ("gradient_entry_" + str(i), covariance.hessian, ((v, 0, 0, 0), 1, 0, 0))
        )
        calls.append(
            ("coefficient_" + str(i), covariance.hessian, ((0, 0, 0, 0), v, 0, 0))
        )
    calls += [
        ("nonpositive_p_" + str(i), covariance.hessian, ((0, 0, 0, 0), v, 0, 0))
        for i, v in enumerate((0, -1))
    ]
    return calls
