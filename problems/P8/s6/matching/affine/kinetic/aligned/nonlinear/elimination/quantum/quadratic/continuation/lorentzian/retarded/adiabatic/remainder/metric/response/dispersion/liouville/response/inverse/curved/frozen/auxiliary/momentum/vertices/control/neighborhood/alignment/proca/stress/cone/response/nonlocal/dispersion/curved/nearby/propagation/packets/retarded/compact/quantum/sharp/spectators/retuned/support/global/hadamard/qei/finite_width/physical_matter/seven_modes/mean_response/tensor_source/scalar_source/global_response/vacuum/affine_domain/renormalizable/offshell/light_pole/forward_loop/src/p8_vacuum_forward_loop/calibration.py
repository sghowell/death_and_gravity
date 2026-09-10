"""Exact actual-parameter error margin, without numerical loop quadrature."""

from fractions import Fraction

import sympy as sp
from p8_polynomial_vacuum import model

from . import angular, leading

DEFAULT_RELATIVE_TOLERANCE = sp.Rational(1, 10**6)


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require a finite exact rational tolerance")
    return sp.Rational(value)


def point(relative_tolerance=DEFAULT_RELATIVE_TOLERANCE):
    eps = rational(relative_tolerance)
    if not 0 < eps < 1:
        raise ValueError("Require a positive relative tolerance strictly below one")
    p = model.data()["actual_parameters"]
    M = p["heavy_mass_squared"]
    g = p["cubic_coupling_squared"]
    c = g / M**2
    tree = 4 * p["lambda"]
    angular_error = angular.data()["actual_angular_b2_error_upper"]
    leading_error = leading.data()["actual_angular_zero_b2_upper"]
    error = angular_error + leading_error
    return {
        "relative_tolerance": eps,
        "actual_M": M,
        "actual_cubic_squared": g,
        "actual_small_c": c,
        "actual_tree_b2": tree,
        "angular_b2_error_upper": angular_error,
        "angular_zero_b2_error_upper": leading_error,
        "total_one_loop_b2_error_upper": error,
        "relative_one_loop_error_upper": error / tree,
        "requested_tolerance_proved": bool(error < eps * tree),
        "tree_plus_one_loop_b2_lower": tree - error,
        "scope": "Full one-loop coefficient in the declared fixed subtraction scheme only; higher-loop errors, high-energy contour, finite gravity and common bounce matching are not bounded.",
    }


def data():
    p = model.data()["actual_parameters"]
    M = p["heavy_mass_squared"]
    D = p["D"]
    g = p["cubic_coupling_squared"]
    c = g / M**2
    l = p["lambda"]
    quartic = p["bare_polynomial_quartic"]
    d = p["positive_completed_square_quartic_margin"]
    q = point()
    return {
        "actual_parameters": p,
        "actual_margin": q,
        "checks": {
            "actual_cubic_fixed_interaction_relation": g - 2 * l * D**3,
            "actual_full_tree_forward_coefficient": 2 * g / D**3 - 4 * l,
            "actual_small_effective_contact_ratio": sp.factor(
                d / c - (4 + 4 / D - 8 / D**2)
            ),
            "actual_angular_to_tree_ratio": sp.factor(
                q["angular_b2_error_upper"] / (4 * l) - 30 * g * D**3 / M**3
            ),
        },
        "bounds": {
            "small_effective_contact_below_five_c": 0 < d < 5 * c,
            "uniform_C_bound_below_five_c_M": quartic + g / (M - 3) < 5 * c * M,
            "B_bound_eleven_c_sufficient": 5 + 3 * M / (M - 3) < 11,
            "actual_complete_one_loop_error_below_one_millionth_tree": q[
                "requested_tolerance_proved"
            ],
            "tree_plus_complete_one_loop_coefficient_strictly_positive": q[
                "tree_plus_one_loop_b2_lower"
            ]
            > 0,
            "leading_error_negligible_but_not_dropped": 0
            < q["angular_zero_b2_error_upper"] / (4 * l)
            < sp.Rational(1, 10**200),
            "angular_error_strictly_below_thirty_cubic_squared_relative": q[
                "angular_b2_error_upper"
            ]
            / (4 * l)
            < 30 * g,
            "small_couplings_alone_do_not_prove_relative_matching": quartic
            * quartic
            / (4 * l)
            > 10**180,
        },
    }


def bad_cases():
    invalid = (
        True,
        False,
        1.0,
        sp.Float(1),
        "1",
        sp.I,
        sp.oo,
        -sp.oo,
        sp.zoo,
        sp.nan,
        None,
        sp.Symbol("epsilon"),
    )
    rows = [
        ("inexact_or_nonreal_tolerance_" + str(i), point, (v,))
        for i, v in enumerate(invalid)
    ]
    rows += [
        ("outside_relative_tolerance_domain_" + str(i), point, (v,))
        for i, v in enumerate((0, -1, 1, 2))
    ]
    return rows
