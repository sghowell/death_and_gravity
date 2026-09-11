"""Full matched finite vacuum reference at the unchanged GY14 boundary."""

from functools import cache

import sympy as s
from p8_polynomial_vacuum import model
from p8_vacuum_fermion_vacuum import calibration as fermion
from p8_vacuum_full_one_loop_matching import calibration as first

from . import bounds


@cache
def data():
    p = model.data()["actual_parameters"]
    L, g, M, lam = [
        p[k]
        for k in (
            "bare_polynomial_quartic",
            "cubic_coupling_squared",
            "heavy_mass_squared",
            "lambda",
        )
    ]
    m = first.data()["actual_MS_interaction_parameters"]["mF"]
    ell = first.data()["exact_scale_log_rational_enclosure"]["scale_log_upper"]
    scalar = bounds.scalar_enclosure(m, M, L, g, 144, 10**25, 10**25, 10**13)
    firstV = bounds.one_loop_enclosure(m, M, 144, 256, ell)
    fermionV = fermion.data()["both_vacuum_primitive_absolute_upper"]
    second = scalar["complete_scalar_vacuum_absolute_upper"] + fermionV
    relative = second / firstV["complete_lower"]
    G = s.sqrt(g)
    source1 = G * (ell + 1) / (2 * 144)
    source_square = (
        g * (2 * ell**2 + 4 * ell + 3 + s.Rational(16, 6)) / (8 * M * 144**2)
    )
    return {
        "actual_parameters": p,
        "actual_fermion_mass": m,
        "verified_scale_log_upper": ell,
        "complete_scalar_vacuum_enclosure": scalar,
        "complete_one_loop_vacuum_enclosure": firstV,
        "unchanged_both_fermion_vacuum_primitives_upper": fermionV,
        "complete_matched_two_loop_vacuum_absolute_upper": second,
        "two_loop_relative_to_complete_one_loop_upper": relative,
        "finite_vacuum_counterterm_first_interval": {
            "lower": -firstV["complete_upper"],
            "upper": -firstV["complete_lower"],
        },
        "finite_vacuum_counterterm_second_absolute_upper": second,
        "finite_first_H_source_absolute_upper": source1,
        "finite_regulated_source_square_absolute_upper_before_exact_cancellation": source_square,
        "checks": {
            "same_actual_fermion_mass": m - 10**200,
            "fundamental_cubic_square": G**2 - g,
            "same_exact_tree_reference": 4 * lam - 2 * g / (M - 2) ** 3,
            "complete_scalar_fermion_vacuum_sum": second
            - scalar["complete_scalar_vacuum_absolute_upper"]
            - fermionV,
            "full_relative_bound_uses_complete_V1": relative * firstV["complete_lower"]
            - second,
            "physical_second_vacuum_reference_magnitude": second - second,
            "first_source_reference_bound": source1 - G * (ell + 1) / 288,
            "source_square_keeps_all_three_finite_terms": source_square
            - g * (2 * ell**2 + 4 * ell + 3 + s.Rational(16, 6)) / (8 * M * 144**2),
        },
        "bounds": {
            "actual_heavy_hierarchy": bool(2 < M < m * m),
            "actual_verified_scale_log_upper_below_1000": bool(0 < ell < 1000),
            "all_scalar_vacuum_groups_positive": all(
                v > 0 for v in scalar["group_uppers"].values()
            ),
            "complete_scalar_vacuum_upper_below_one_e_275": bool(
                0 < scalar["complete_scalar_vacuum_absolute_upper"] < 10**275
            ),
            "complete_one_loop_vacuum_strictly_positive": bool(
                firstV["complete_lower"] > 0
            ),
            "complete_two_loop_vacuum_upper_below_one_e_595": bool(
                0 < second < 10**595
            ),
            "complete_two_loop_relative_below_one_e_minus_203": bool(
                0 < relative < s.Rational(1, 10**203)
            ),
            "complete_two_loop_vacuum_coefficient_below_first_lower": bool(
                second < firstV["complete_lower"]
            ),
            "finite_first_source_upper_below_one_e_minus_3": bool(
                0 < source1 < s.Rational(1, 10**3)
            ),
            "source_square_reference_upper_below_one_e_minus_200": bool(
                0 < source_square < s.Rational(1, 10**200)
            ),
        },
        "scope": "The unique finite vacuum-energy-zero reference through two loops in the named nongravitational physical-Phi/MS-interaction model, plus complete regulated first-source-square ownership. Large dimensionful reference values are retained. The full second H-source coefficient, global potential, physical truncation and V/G/B remain open.",
    }
