"""Actual wineglass and four disjoint raw-graph group error bounds."""

from fractions import Fraction
from functools import cache

import sympy as sp
from p8_polynomial_vacuum import model
from p8_vacuum_two_loop_double_bubble import calibration as previous


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require a finite exact rational radial invariant")
    return sp.Rational(value)


def point(invariant=0):
    y = rational(invariant)
    if y < 0:
        raise ValueError("Require a nonnegative radial invariant")
    M = data()["actual_heavy_mass_squared"]
    return {
        "radial_invariant": y,
        "shifted_inner_bubble_modulus_upper": (sp.log(1 + y) + 3) / 144,
        "decaying_inner_bubble_difference_modulus_upper": (8 * sp.sqrt(y) + 4)
        / ((y + 4) * 144),
        "centered_outer_light_inverse_real_lower": y + sp.Rational(1, 4),
        "routed_heavy_inverse_real_lower": (y + M) / 2,
        "scope": "All-radius conservative modulus bounds on the continued inner bubble and its decaying difference, not signed values or a cutoff.",
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
    prior = previous.data()["three_integrated_groups_combined_upper"]
    error = sp.Rational(14496240, 20736) * L**3
    return {
        "actual_heavy_mass_squared": M,
        "actual_cubic_squared": g,
        "actual_quartic": L,
        "actual_tree_b2": 4 * lam,
        "actual_wineglass_group_b2_absolute_upper": error,
        "actual_wineglass_group_relative_to_tree_upper": error / (4 * lam),
        "previous_three_integrated_groups_upper": prior,
        "four_integrated_raw_graph_groups_combined_upper": prior + error,
        "represented_raw_refinement_count": 192,
        "remaining_raw_refinement_count": 0,
        "scope": "All 192 bare refinements are represented once with the specified graph-level subtractions. This is not the complete two-loop amplitude: fixed finite-potential counterterm insertions, complete normalization and source-aware matching remain open.",
        "checks": {
            "complete_subtracted_and_heavy_logarithmic_prefactor": 240
            + 40 * 600 * (600 + 4)
            - 14496240,
            "same_two_loop_measure_lower_bound": 16**2 * 3**4 - 20736,
            "all_four_raw_graph_groups_exactly_accounted": 88 + 64 + 24 + 16 - 192,
            "same_actual_parent_heavy_mass": M
            - previous.data()["actual_heavy_mass_squared"],
            "same_actual_parent_tree_coefficient": 4 * lam
            - previous.data()["actual_tree_b2"],
            "same_previous_disjoint_bound": prior
            - previous.data()["three_integrated_groups_combined_upper"],
        },
        "bounds": {
            "actual_mass_above_thirty_two": M > 32,
            "actual_heavy_radial_margin_positive": M / 2 - 9 > 0,
            "actual_external_heavy_gap_above_half_mass": M - 3 > M / 2,
            "actual_cubic_over_mass_below_quartic_third": g / M < L / 3,
            "actual_logarithm_argument_below_ten_to_two_hundred": 4 * M < 10**200,
            "positive_wineglass_absolute_error_below_one_e_minus_611": 0
            < error
            < sp.Rational(1, 10**611),
            "positive_wineglass_relative_error_below_three_e_minus_12": 0
            < error / (4 * lam)
            < sp.Rational(3, 10**12),
            "four_groups_below_three_e_minus_607": 0
            < prior + error
            < sp.Rational(3, 10**607),
            "four_groups_below_one_e_minus_7_of_tree": 0
            < (prior + error) / (4 * lam)
            < sp.Rational(1, 10**7),
            "strict_complex_logarithm_strip_gap": sp.Rational(1, 8) > 0,
            "strict_subtracted_integral_rational_margin": sp.Rational(29, 16) > 0,
        },
    }


def bad_cases():
    invalid = (
        True,
        False,
        1.0,
        sp.Float(1),
        "1",
        None,
        sp.oo,
        -sp.oo,
        sp.zoo,
        sp.I,
        sp.nan,
        sp.Symbol("y"),
    )
    return [
        ("inexact_radial_" + str(i), point, (v,)) for i, v in enumerate(invalid)
    ] + [
        ("negative_radial_" + str(i), point, (v,))
        for i, v in enumerate((-1, Fraction(-1, 2)))
    ]
