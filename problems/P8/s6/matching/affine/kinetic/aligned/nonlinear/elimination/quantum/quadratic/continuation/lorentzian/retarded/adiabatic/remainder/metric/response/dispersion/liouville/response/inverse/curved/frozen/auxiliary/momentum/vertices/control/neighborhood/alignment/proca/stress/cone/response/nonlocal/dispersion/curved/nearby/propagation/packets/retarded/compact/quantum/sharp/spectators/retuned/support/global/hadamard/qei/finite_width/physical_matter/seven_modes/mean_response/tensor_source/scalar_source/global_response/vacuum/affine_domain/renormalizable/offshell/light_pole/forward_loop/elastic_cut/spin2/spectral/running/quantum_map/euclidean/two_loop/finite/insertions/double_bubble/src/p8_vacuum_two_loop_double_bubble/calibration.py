"""Exact actual-coupling comparison for the disjoint 24-refinement group."""

from fractions import Fraction
from functools import cache

import sympy as sp
from p8_polynomial_vacuum import model
from p8_vacuum_two_loop_insertions import calibration as previous

from . import bubble


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require a finite exact rational channel invariant")
    return sp.Rational(value)


def point(invariant=0):
    z = rational(invariant)
    if abs(z) > 3:
        raise ValueError("Require a real diagnostic channel invariant in [-3,3]")
    d = data()
    M, g, L = (
        d["actual_heavy_mass_squared"],
        d["actual_cubic_squared"],
        d["actual_quartic"],
    )
    return {
        "channel_invariant": z,
        "actual_external_pair_vertex": -L + g / (M - z),
        "external_pair_vertex_modulus_upper": 2 * L,
        "finite_bubble_modulus_upper": sp.Rational(1, 72),
        "each_summed_heavy_triangle_modulus_upper": 8 * g * 600 / (144 * M),
        "scope": "Exact diagnostic value for the external vertex and conservative modulus bounds for the integrated finite factors; not a signed two-loop amplitude.",
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
    error = sp.Rational(19224, 144**2) * L**3
    prior = previous.data()["two_integrated_groups_combined_upper"]
    ell, a, c = sp.symbols(
        "logarithm_upper positive_quartic positive_cubic_over_mass", positive=True
    )
    raw_prefactor = sp.Rational(3, 4) * (
        (2 * a) ** 3 * 2**2 + 2 * (2 * a) ** 2 * 2 * 8 * c * ell
    )
    return {
        "actual_heavy_mass_squared": M,
        "actual_cubic_squared": g,
        "actual_quartic": L,
        "actual_tree_b2": 4 * lam,
        "actual_double_bubble_group_b2_absolute_upper": error,
        "actual_double_bubble_group_relative_to_tree_upper": error / (4 * lam),
        "previous_two_integrated_groups_upper": prior,
        "three_integrated_groups_combined_upper": prior + error,
        "represented_raw_refinement_count": 176,
        "remaining_raw_wineglass_refinement_count": 16,
        "scope": "Three disjoint raw-graph groups with their specified subtractions, not a complete two-loop amplitude. The 16 wineglass refinements, additional finite-potential counterterm insertions and complete two-loop pole/LSZ accounting remain open.",
        "checks": {
            "all_channel_vertex_loop_and_triangle_prefactors": sp.expand(
                raw_prefactor.subs(c, a / 3) - (24 + 32 * ell) * a**3
            ),
            "fixed_elementary_logarithm_majorant_prefactor": 24 + 32 * 600 - 19224,
            "two_loop_measure_lower_uses_pi_above_three": 16**2 * 3**4 - 144**2,
            "three_raw_groups_are_disjoint_and_leave_sixteen": 88 + 64 + 24 + 16 - 192,
            "same_previous_two_group_bound": prior
            - previous.data()["two_integrated_groups_combined_upper"],
            "logarithm_four_bound_has_positive_exponential_witness": bubble.data()[
                "positive_exp_two_partial_sum"
            ]
            - 5,
        },
        "bounds": {
            "actual_mass_above_thirty_two": M > 32,
            "actual_external_heavy_gap_above_half_mass": M - 3 > M / 2,
            "actual_internal_heavy_radial_margin_positive": M / 2 - 9 > 0,
            "actual_triangle_mass_ratio_below_two": M**2
            < 2 * (M - sp.Rational(1, 4)) ** 2,
            "actual_cubic_over_mass_below_quartic_third": g / M < L / 3,
            "actual_logarithm_argument_below_ten_to_two_hundred": 4 * M < 10**200,
            "elementary_log_four_below_two_witness": bubble.data()[
                "positive_exp_two_partial_sum"
            ]
            > 4,
            "positive_group_absolute_error_below_two_e_minus_614": 0
            < error
            < sp.Rational(2, 10**614),
            "positive_group_relative_error_below_three_e_minus_15": 0
            < error / (4 * lam)
            < sp.Rational(3, 10**15),
            "three_groups_below_three_e_minus_607": 0
            < prior + error
            < sp.Rational(3, 10**607),
            "three_groups_below_one_e_minus_7_of_tree": 0
            < (prior + error) / (4 * lam)
            < sp.Rational(1, 10**7),
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
        sp.Symbol("z"),
    )
    return [
        ("inexact_channel_" + str(i), point, (v,)) for i, v in enumerate(invalid)
    ] + [
        ("outside_channel_" + str(i), point, (v,))
        for i, v in enumerate((-4, 4, Fraction(7, 2)))
    ]
