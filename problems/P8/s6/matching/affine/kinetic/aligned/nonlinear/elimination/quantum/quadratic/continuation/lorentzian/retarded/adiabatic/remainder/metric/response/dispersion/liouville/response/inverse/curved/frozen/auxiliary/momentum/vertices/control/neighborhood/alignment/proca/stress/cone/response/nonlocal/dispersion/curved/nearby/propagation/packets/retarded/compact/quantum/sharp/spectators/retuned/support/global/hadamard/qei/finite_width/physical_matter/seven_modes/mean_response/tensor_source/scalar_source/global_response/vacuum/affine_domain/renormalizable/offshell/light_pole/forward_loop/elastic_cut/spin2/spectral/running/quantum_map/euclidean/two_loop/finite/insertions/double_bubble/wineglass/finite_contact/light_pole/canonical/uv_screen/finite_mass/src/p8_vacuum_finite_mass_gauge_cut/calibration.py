"""Exact finite-mass relative bound at the same prospective GY14 boundary."""

from functools import cache

import sympy as sp
from p8_vacuum_gauge_yukawa_screen import calibration as parent

from . import cuts


@cache
def data():
    inherited = parent.data()
    m = inherited["new_common_Dirac_mass"]
    a = inherited["prospective_gauge_squared"]
    Y = inherited["prospective_Yukawa_squared"]
    z = a * Y / (16 * sp.pi**2 * m * m)
    interior = cuts.point(mass=m)
    eps = interior["per_box_error_over_z_upper"]
    relative = interior["relative_boundary_difference_error_upper"]
    strict_upper = sp.Rational(5, 10**194)
    checks = {
        "same_fermion_mass": m - 10**200,
        "same_positive_transition_scale": sp.simplify(
            z + 3 * inherited["leading_threshold_coefficient_at_prospective_boundary"]
        ),
        "actual_Neumann_ratio": interior["Neumann_ratio_upper"]
        - sp.Rational(12, 10**200),
        "actual_per_box_relative_error": eps - sp.Rational(1, 10**194),
        "actual_cut_relative_bound_dictionary": relative
        - sp.Rational(9, 2) * (eps + eps * eps / 8),
        "strict_cut_relative_upper_dictionary": strict_upper - 5 * eps,
        "dimensionless_leading_interior_jump": interior[
            "leading_boundary_difference_over_dA_z2_div_9pi"
        ]
        + 4,
    }
    bounds = {
        "same_positive_finite_fermion_mass": m >= 24,
        "actual_strict_Neumann_disc": 0 < 12 / m < sp.Rational(1, 2),
        "gauge_and_Yukawa_parameters_positive": bool(a > 0 and Y > 0),
        "positive_nonzero_transition_scale": z > 0,
        "per_box_error_below_half": 0 < eps < sp.Rational(1, 2),
        "finite_mass_relative_cut_error_below_five_e_minus_194": 0
        < relative
        < strict_upper,
        "finite_mass_relative_cut_error_below_one": relative < 1,
        "strict_nonzero_interior_two_gauge_boundary_mismatch": interior[
            "strict_nonzero_boundary_difference_at_this_order"
        ],
        "opposite_interior_channel_also_nonzero": cuts.point(sp.Rational(5, 2), m)[
            "strict_nonzero_boundary_difference_at_this_order"
        ],
        "crossing_center_not_falsely_certified_nonzero": not cuts.point(2, m)[
            "strict_nonzero_boundary_difference_at_this_order"
        ],
        "small_mass_bound_can_be_inconclusive": not cuts.point(mass=24)[
            "strict_nonzero_boundary_difference_at_this_order"
        ],
        "near_center_bound_can_be_inconclusive": not cuts.point(
            2 + sp.Rational(1, 10**200), m
        )["strict_nonzero_boundary_difference_at_this_order"],
    }
    return {
        "candidate": "GY14-unbroken, exactly the prospective S6.128 boundary",
        "fermion_mass": m,
        "box_scale_z": z,
        "per_box_error_over_z_upper": eps,
        "first_two_gauge_cut_relative_error_upper_at_interior": relative,
        "strict_simple_relative_upper": strict_upper,
        "interior": interior,
        "center": cuts.point(2, m),
        "small_mass_control": cuts.point(mass=24),
        "near_center_control": cuts.point(2 + sp.Rational(1, 10**200), m),
        "checks": checks,
        "bounds": bounds,
        "scope": "This resolves the one-loop fermion-box momentum remainder needed to retain the first two-gauge cut in S6.128. It does not compute higher-loop matching, a new-model complete b2 budget, confinement, finite-gravity contours or a common bounce parent.",
    }
