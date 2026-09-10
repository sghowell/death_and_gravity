"""Actual grouped insertion error, retaining the fixed one-loop comparison."""

from fractions import Fraction
from functools import cache

import sympy as sp
from p8_polynomial_vacuum import model
from p8_vacuum_euclidean_control import calibration as euclidean
from p8_vacuum_forward_loop import calibration as one_loop
from p8_vacuum_two_loop_finite import calibration as finite

from . import radial


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require a finite exact rational radial invariant")
    return sp.Rational(value)


def point(invariant=0):
    y = rational(invariant)
    if y < 0:
        raise ValueError("Require a nonnegative centered Euclidean radial invariant")
    d = data()
    M, g = d["actual_heavy_mass_squared"], d["actual_cubic_squared"]
    logarithm = 1 / M if y == 0 else sp.log(1 + y / M) / y
    return {
        "radial_invariant": y,
        "decaying_remainder_logarithmic_upper": 2 * g * logarithm / 144,
        "decaying_remainder_uniform_upper": g / (144 * M),
        "combined_light_denominator_real_lower": y + sp.Rational(1, 4),
        "scope": "Upper bounds for the modulus of the decaying kernel remainder along every routed complex shift; these are not signed self-energy values.",
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
    alpha = euclidean.data()["actual_uniform_relative_self_energy_upper"]
    E1 = one_loop.point()["total_one_loop_b2_error_upper"]
    constant = 2 * alpha * E1
    remainder = 54 * L**2 * g * 600 / (144**2 * (M - sp.Rational(1, 4)))
    error = constant + remainder
    previous = finite.data()["actual_finite_subsector_b2_absolute_upper"]
    r = radial.data()
    Mvar = r["M"]
    return {
        "actual_heavy_mass_squared": M,
        "actual_cubic_squared": g,
        "actual_quartic": L,
        "actual_tree_b2": 4 * lam,
        "same_fixed_asymptotic_multiplier_upper": alpha,
        "same_complete_renormalized_one_loop_b2_upper": E1,
        "constant_insertion_group_b2_upper": constant,
        "decaying_insertion_group_b2_upper": remainder,
        "complete_grouped_insertion_b2_absolute_upper": error,
        "grouped_insertion_relative_to_tree_upper": error / (4 * lam),
        "previous_UV_finite_group_b2_upper": previous,
        "two_integrated_groups_combined_upper": previous + error,
        "raw_refinements_represented_in_both_groups": 152,
        "remaining_subtraction_dependent_raw_refinements": 40,
        "scope": "Both integrated groups are still only part of the full two-loop four-point amplitude. The 64 tadpole refinements are grouped with their fixed inner OS and linear inherited outer counterterms; 40 other raw refinements and complete two-loop pole/LSZ normalization remain open.",
        "checks": {
            "same_once_fixed_kinetic_majorant": alpha - g / (288 * M),
            "same_parent_one_loop_majorant": E1
            - one_loop.point()["total_one_loop_b2_error_upper"],
            "constant_group_uses_two_lines": constant - 2 * alpha * E1,
            "all_raw_refinement_groups_accounted": 88 + 64 + 40 - 192,
            "same_positive_radial_mass_gap": r["actual_gap_one_quarter_integral"]
            - sp.log(4 * Mvar) / (Mvar - sp.Rational(1, 4)),
        },
        "bounds": {
            "actual_mass_above_thirty_two": M > 32,
            "heavy_denominator_gap_above_half_mass": M - 9 > M / 2,
            "routed_kernel_strip_denominator_gap_above_half": sp.Rational(1, 2) - 2 / M
            > 0,
            "actual_cubic_over_mass_below_quartic_third": g / M < L / 3,
            "actual_logarithm_argument_below_ten_to_two_hundred": 4 * M < 10**200,
            "grouped_absolute_error_below_one_e_minus_614": 0
            < error
            < sp.Rational(1, 10**614),
            "grouped_relative_error_below_two_e_minus_15": 0
            < error / (4 * lam)
            < sp.Rational(2, 10**15),
            "constant_piece_retained_below_one_e_minus_813": 0
            < constant
            < sp.Rational(1, 10**813),
            "two_integrated_groups_below_three_e_minus_607": 0
            < previous + error
            < sp.Rational(3, 10**607),
            "two_integrated_groups_below_one_e_minus_7_of_tree": 0
            < (previous + error) / (4 * lam)
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
        sp.Symbol("y"),
    )
    rows = [("inexact_radial_" + str(i), point, (v,)) for i, v in enumerate(invalid)]
    rows += [
        ("negative_radial_" + str(i), point, (v,))
        for i, v in enumerate((-1, Fraction(-1, 2)))
    ]
    return rows
