"""Actual canonical b2 bounds at tree, one loop and two loops."""

from functools import cache

import sympy as sp
from p8_vacuum_forward_loop import calibration as one_loop
from p8_vacuum_two_loop_finite_contact import calibration as previous
from p8_vacuum_two_loop_light_pole import calibration as pole


def point(order=2):
    if type(order) is not int:
        raise TypeError("Require a native integer perturbative order")
    if not 0 <= order <= 2:
        raise ValueError("Only orders zero, one and two are certified here")
    d = data()
    error = (
        sp.S.Zero,
        d["actual_one_loop_b2_absolute_upper"],
        d["actual_one_plus_two_loop_b2_absolute_upper"],
    )[order]
    tree = d["actual_tree_b2"]
    return {
        "highest_loop_order": order,
        "canonical_b2_absolute_error_upper": error,
        "canonical_b2_lower": tree - error,
        "canonical_b2_upper": tree + error,
        "canonical_b2_relative_error_upper": error / tree,
        "strictly_positive_canonical_b2_at_this_order": bool(tree - error > 0),
        "scope": "The gravitationally decoupled scalar canonical Phi elastic forward coefficient in the fixed scheme through the stated loop order. This is not the full all-orders observable, a two-loop derivative-coordinate/source result or a UV-admissibility verdict.",
    }


@cache
def data():
    d = previous.data()
    tree = d["actual_tree_b2"]
    E1 = one_loop.point()["total_one_loop_b2_error_upper"]
    E2 = d["raw_graph_and_finite_contact_two_loop_upper"]
    pole_bound = pole.data()["actual_one_plus_two_loop_quadratic_coefficient_upper"]
    total = E1 + E2
    return {
        "actual_tree_b2": tree,
        "actual_one_loop_b2_absolute_upper": E1,
        "actual_complete_two_loop_b2_absolute_upper": E2,
        "actual_one_plus_two_loop_b2_absolute_upper": total,
        "actual_complete_two_loop_relative_error_upper": E2 / tree,
        "actual_one_plus_two_loop_relative_error_upper": total / tree,
        "actual_two_loop_canonical_b2_lower": tree - total,
        "actual_first_sheet_unit_disc_light_inverse_error_upper": pole_bound,
        "scope": "The previous raw-graph and finite-contact sum now bounds the full canonical two-loop b2 only after the independent marked-core, local-reference, cubic-product and external-field ledger is established. The exact rational estimates do not determine the correction's sign, all higher-loop errors, all-energy cuts/contours or original P8 closure.",
        "checks": {
            "same_once_fixed_actual_tree_coefficient": tree
            - previous.data()["actual_tree_b2"],
            "same_actual_complete_group_bound": E2
            - previous.data()["raw_graph_and_finite_contact_two_loop_upper"],
            "same_inherited_known_one_plus_two_loop_sum": total
            - d["known_one_and_two_loop_contributions_upper"],
            "same_positive_tree_minus_error_lower": tree
            - total
            - d["positive_tree_minus_known_corrections_lower"],
        },
        "bounds": {
            "actual_tree_b2_positive": tree > 0,
            "actual_one_loop_error_positive": E1 > 0,
            "actual_complete_two_loop_error_positive": E2 > 0,
            "actual_complete_two_loop_error_below_three_e_minus_607": E2
            < sp.Rational(3, 10**607),
            "actual_complete_two_loop_relative_error_below_one_e_minus_7": E2 / tree
            < sp.Rational(1, 10**7),
            "actual_one_plus_two_loop_relative_error_below_one_e_minus_6": total / tree
            < sp.Rational(1, 10**6),
            "actual_two_loop_canonical_b2_strictly_positive": tree - total > 0,
            "actual_light_pole_normalization_control_below_one_e_minus_18": 0
            < pole_bound
            < sp.Rational(1, 10**18),
        },
    }


def bad_cases():
    bad = (
        True,
        False,
        1.0,
        sp.Integer(1),
        sp.Float(1),
        "1",
        None,
        sp.oo,
        -sp.oo,
        sp.I,
        sp.nan,
        sp.Symbol("order"),
    )
    return [("inexact_order_" + str(i), point, (v,)) for i, v in enumerate(bad)] + [
        ("unsupported_order_" + str(i), point, (v,)) for i, v in enumerate((-1, 3, 100))
    ]
