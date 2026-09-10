"""Actual coupling-weighted finite-subsector enclosures without floating logs."""

from fractions import Fraction
from functools import cache
from math import factorial

import sympy as sp
from p8_polynomial_vacuum import model
from p8_vacuum_two_loop_denom import graphs, subgraphs

from . import ordered, parametric, sectors, selection


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require a finite exact rational logarithm upper bound")
    return sp.Rational(value)


def polynomial(kind, choices, log_upper=600):
    selection.require_finite(kind, choices)
    ell = rational(log_upper)
    if ell < 0:
        raise ValueError("Require a nonnegative logarithm upper bound")
    coefficients = sectors.data(kind, choices)["box_log_polynomial_coefficients"]
    value = sum(
        sp.Rational(c.numerator, c.denominator) * ell**j
        for j, c in enumerate(coefficients)
    )
    return {
        "kind": kind,
        "choices": choices,
        "log_upper": ell,
        "scaled_anisotropic_box_integral_upper": value,
        "scope": "A conditional polynomial upper bound when log B is no larger than the supplied value. Only the actual-parameter calibration below proves log B below 600.",
    }


@cache
def data():
    p = model.data()["actual_parameters"]
    M, g, L, lam = (
        p[k]
        for k in (
            "heavy_mass_squared",
            "cubic_coupling_squared",
            "bare_polynomial_quartic",
            "lambda",
        )
    )
    coefficients = [sp.Integer(0)] * 3
    weights = graphs.data()["family_Wick_weights"]
    for row in sectors.summaries():
        h = row["heavy_edge_count"]
        pref = parametric.data()["per_heavy_edge_count"][h][
            "rational_coupling_and_loop_prefactor"
        ]
        for j, c in enumerate(row["box_log_polynomial_coefficients"]):
            coefficients[j] += (
                weights[row["kind"]] * pref * sp.Rational(c.numerator, c.denominator)
            )
    C = sum(c * 600**j for j, c in enumerate(coefficients))
    error = sp.factor(C * L**3)
    D = sp.Symbol("positive_D", positive=True)
    exp_partial = sum(sp.Rational(3) ** j / factorial(j) for j in range(4))
    return {
        "actual_heavy_mass_squared": M,
        "actual_heavy_box_scale": 2 * M,
        "actual_cubic_squared": g,
        "actual_quartic": L,
        "actual_tree_b2": 4 * lam,
        "aggregate_log_polynomial_coefficients": tuple(coefficients),
        "actual_logarithm_strict_upper": sp.Integer(600),
        "aggregate_constant_at_log_upper": C,
        "actual_finite_subsector_b2_absolute_upper": error,
        "actual_finite_subsector_relative_to_tree_upper": error / (4 * lam),
        "elementary_positive_exp_three_partial_sum": exp_partial,
        "scope": "An absolute error enclosure for the sum of the 88 UV-finite bare two-loop four-point refinements and all external assignments. This is not the signed correction or the complete two-loop coefficient.",
        "checks": {
            "actual_coupling_ratio_positive_gap": sp.factor(
                (3 * D - 2) * (D + 2) / D**2 - 3 - 4 * (D - 1) / D**2
            ),
            "aggregate_constant_coefficient": coefficients[0]
            - sp.Rational(13688587, 17496),
            "aggregate_linear_log_coefficient": coefficients[1]
            - sp.Rational(515831, 972),
            "aggregate_quadratic_log_coefficient": coefficients[2]
            - sp.Rational(59789, 972),
            "aggregate_at_log_six_hundred": C - sp.Rational(393017383387, 17496),
            "elementary_exp_three_partial_sum": exp_partial - 13,
        },
        "bounds": {
            "actual_mass_above_thirty_two": M > 32,
            "actual_heavy_box_scale_above_sixty_four": 2 * M > 64,
            "actual_heavy_scale_below_ten_to_two_hundred": 2 * M < 10**200,
            "exp_three_positive_partial_sum_above_ten": exp_partial > 10,
            "actual_cubic_over_mass_less_than_quartic_third": g / M < L / 3,
            "finite_subsector_absolute_upper_below_three_e_minus_607": 0
            < error
            < sp.Rational(3, 10**607),
            "finite_subsector_relative_upper_below_one_e_minus_7": 0
            < error / (4 * lam)
            < sp.Rational(1, 10**7),
            "tree_minus_finite_subsector_upper_positive": 4 * lam - error > 0,
        },
    }


def bad_cases():
    base = ("double_bubble", (0, 0, 2))
    rows = [
        (
            "needs_subtractions_" + kind + "_" + "".join(map(str, choices)),
            selection.require_finite,
            (kind, choices),
        )
        for kind, choices in graphs.cases()
        if subgraphs.data(kind, choices)["UV_subgraphs"]
    ]
    rows += [
        ("inexact_log_" + str(i), polynomial, base + (v,))
        for i, v in enumerate(
            (True, False, 1.0, sp.Float(1), "1", None, sp.oo, sp.I, sp.nan, -1)
        )
    ]
    rows += [
        ("invalid_order_" + str(i), sectors.single, base + (v,))
        for i, v in enumerate(
            (
                [0, 1, 2, 3, 4],
                (0, 1, 2, 3),
                (0, 1, 2, 3, 3),
                (False, 1, 2, 3, 4),
                (0, 1, 2, 3, sp.Integer(4)),
                (0, 1, 2, 3, 5),
            )
        )
    ]
    rows += [
        ("invalid_large_pattern", ordered.large_data, ((0,),)),
        ("boolean_large_pattern", ordered.large_data, ((-2, False),)),
        ("UV_divergent_small_pattern", ordered.small_data, ((-2, 0),)),
        ("boolean_small_pattern", ordered.small_data, ((False,),)),
    ]
    return rows
