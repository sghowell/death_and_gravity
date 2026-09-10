"""Exact fixed finite contact and actual inherited amplitude error bounds."""

from fractions import Fraction
from functools import cache

import sympy as sp
from p8_polynomial_vacuum import model
from p8_vacuum_forward_loop import calibration as one_loop
from p8_vacuum_two_loop_wineglass import calibration as previous

from . import radial


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require a finite exact rational radial invariant")
    return sp.Rational(value)


def point(invariant=0):
    y = rational(invariant)
    if y < 0:
        raise ValueError("Require a nonnegative radial invariant")
    d = data()
    M, g, L = (
        d["actual_heavy_mass_squared"],
        d["actual_cubic_squared"],
        d["actual_quartic"],
    )
    A = (L - g / M) / 2
    F = A - g / (y + M)
    return {
        "radial_invariant": y,
        "actual_positive_F": F,
        "actual_positive_A": A,
        "negative_finite_contact_kernel": F**2 - A**2,
        "finite_contact_density_without_loop_measure": 6
        * y
        * (F**2 - A**2)
        / (y + 1) ** 2,
        "scope": "Exact all-radius integrand data for the once-fixed finite potential contact. The signed integrated contact is fixed by its anchored integral, not fitted to the forward coefficient.",
    }


@cache
def data():
    p = model.data()["actual_parameters"]
    M, g, L, lam, margin = (
        p[k]
        for k in (
            "heavy_mass_squared",
            "cubic_coupling_squared",
            "bare_polynomial_quartic",
            "lambda",
            "positive_completed_square_quartic_margin",
        )
    )
    r = radial.data()
    J1, J2 = (r[k].subs(r["M"], M) for k in ("J1", "J2"))
    C0 = -L + g / M
    contact = 6 * (C0 * g * J1 + g**2 * J2) / (16 * sp.pi**2)
    contact_upper = sp.Rational(50, 3) * L**2
    derivative_upper = sp.Rational(4812, 144) * L
    error = contact_upper * derivative_upper
    prior = previous.data()["four_integrated_raw_graph_groups_combined_upper"]
    E1 = one_loop.point()["total_one_loop_b2_error_upper"]
    known_two = prior + error
    return {
        "actual_heavy_mass_squared": M,
        "actual_cubic_squared": g,
        "actual_quartic": L,
        "actual_tree_b2": 4 * lam,
        "actual_completed_square_quartic_margin": margin,
        "actual_fixed_finite_potential_contact": contact,
        "actual_contact_absolute_upper": contact_upper,
        "actual_one_loop_quartic_derivative_b2_upper": derivative_upper,
        "actual_finite_contact_insertion_b2_absolute_upper": error,
        "actual_finite_contact_insertion_relative_to_tree_upper": error / (4 * lam),
        "previous_four_raw_graph_groups_upper": prior,
        "raw_graph_and_finite_contact_two_loop_upper": known_two,
        "known_one_and_two_loop_contributions_upper": E1 + known_two,
        "positive_tree_minus_known_corrections_lower": 4 * lam - E1 - known_two,
        "scope": "The once-fixed finite potential contact and its one-loop insertion are included with all four raw-graph groups. These quantitative fixed-order contributions are not by themselves a completed two-loop pole/residue/LSZ or source-aware matching calculation, an all-orders truncation bound or original P8 closure.",
        "checks": {
            "fixed_contact_from_same_zero_external_amplitude": contact
            - 6 * (C0 * g * J1 + g**2 * J2) / (16 * sp.pi**2),
            "same_constant_potential_high_momentum_vertex": C0 + L - g / M,
            "exact_rational_contact_majorant_prefactor": sp.Rational(4 * 600, 144)
            - sp.Rational(50, 3),
            "exact_rational_derivative_majorant_prefactor": sp.Rational(
                12 + 8 * 600, 144
            )
            - sp.Rational(4812, 144),
            "exact_product_insertion_majorant": error - sp.Rational(10025, 18) * L**3,
            "same_parent_raw_graph_bound": prior
            - previous.data()["four_integrated_raw_graph_groups_combined_upper"],
        },
        "bounds": {
            "actual_mass_above_thirty_two": M > 32,
            "actual_F_zero_strictly_positive": (L - 3 * g / M) / 2 > 0,
            "actual_C_zero_strictly_negative": C0 < 0,
            "actual_cubic_over_mass_below_quartic_third": g / M < L / 3,
            "actual_logarithm_argument_below_ten_to_two_hundred": 4 * M < 10**200,
            "actual_mass_ratio_majorant_below_two": M**2 < 2 * (M - 1) ** 2,
            "contact_upper_below_one_e_minus_408": 0
            < contact_upper
            < sp.Rational(1, 10**408),
            "finite_contact_upper_below_one_e_minus_6_of_tree_quartic_margin": 0
            < contact_upper / margin
            < sp.Rational(1, 10**6),
            "insertion_error_below_seven_e_minus_612": 0
            < error
            < sp.Rational(7, 10**612),
            "insertion_relative_error_below_two_e_minus_12": 0
            < error / (4 * lam)
            < sp.Rational(2, 10**12),
            "known_two_loop_groups_below_three_e_minus_607": 0
            < known_two
            < sp.Rational(3, 10**607),
            "known_two_loop_groups_below_one_e_minus_7_of_tree": 0
            < known_two / (4 * lam)
            < sp.Rational(1, 10**7),
            "known_one_plus_two_loop_groups_below_one_e_minus_6_of_tree": 0
            < (E1 + known_two) / (4 * lam)
            < sp.Rational(1, 10**6),
            "tree_minus_known_correction_bound_strictly_positive": 4 * lam
            - E1
            - known_two
            > 0,
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
