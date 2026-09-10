"""Exact fixed-order spin-two slope and negative-transfer limit calibrations."""

from fractions import Fraction

import sympy as sp
from p8_polynomial_vacuum import model

from . import triangles


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require a finite exact rational")
    return sp.Rational(value)


def point(Planck_mass=10**400, transfer_radius=1):
    P, r = rational(Planck_mass), rational(transfer_radius)
    if P <= 0:
        raise ValueError("Require positive Planck mass")
    if not 0 < r < 4:
        raise ValueError("Require transfer radius strictly between zero and four")
    S = triangles.data()["actual_spin_two_slope_rational_upper"]
    return {
        "Planck_mass_in_light_mass_units": P,
        "transfer_disc_radius": r,
        "spin_two_form_factor_slope_upper": S,
        "subtracted_form_factor_linear_upper_coefficient": S / (1 - r / 4),
        "form_factor_quadratic_remainder_upper_coefficient": S / (4 * (1 - r / 4)),
        "negative_t_vertex_finite_correction_absolute_upper": 2 * S / P**2,
        "negative_t_vertex_limit_error_per_abs_t_upper": S / (2 * P**2 * (1 - r / 4)),
        "scope": "Exact coefficients and bounds at one matter loop and leading inverse-Planck-square order only. Small positive input P is not a claim of perturbative gravitational validity, and no omitted graviton-loop or Regge error is bounded.",
    }


def data():
    d = point()
    p = model.data()["actual_parameters"]
    P = d["Planck_mass_in_light_mass_units"]
    M = p["heavy_mass_squared"]
    g = p["cubic_coupling_squared"]
    return {
        "actual_selected_calibration": d,
        "checks": {
            "actual_finite_vertex_correction_bound": d[
                "negative_t_vertex_finite_correction_absolute_upper"
            ]
            - g / (1296 * M * M * P * P),
            "unit_radius_linear_form_factor_bound": d[
                "subtracted_form_factor_linear_upper_coefficient"
            ]
            - 4 * d["spin_two_form_factor_slope_upper"] / 3,
            "unit_radius_quadratic_form_factor_bound": d[
                "form_factor_quadratic_remainder_upper_coefficient"
            ]
            - d["spin_two_form_factor_slope_upper"] / 3,
            "unit_radius_fixed_t_limit_error_bound": d[
                "negative_t_vertex_limit_error_per_abs_t_upper"
            ]
            - 2 * d["spin_two_form_factor_slope_upper"] / (3 * P * P),
        },
        "bounds": {
            "actual_negative_t_vertex_piece_below_one_e_minus_1205": 0
            < d["negative_t_vertex_finite_correction_absolute_upper"]
            < sp.Rational(1, 10**1205),
            "actual_vertex_piece_to_tree_ratio_below_one_e_minus_605": d[
                "negative_t_vertex_finite_correction_absolute_upper"
            ]
            / (4 * p["lambda"])
            < sp.Rational(1, 10**605),
            "selected_heavy_mass_below_Planck_scale": M < P * P,
            "unit_transfer_disc_within_conservative_analytic_domain": 0
            < d["transfer_disc_radius"]
            < 4,
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
        sp.Symbol("p"),
    )
    rows = []
    for i, v in enumerate(invalid):
        rows += [
            ("inexact_Planck_" + str(i), point, (v,)),
            ("inexact_radius_" + str(i), point, (10**400, v)),
        ]
    rows += [
        ("nonpositive_Planck_" + str(i), point, (v,)) for i, v in enumerate((0, -1))
    ]
    rows += [
        ("outside_transfer_domain_" + str(i), point, (10**400, v))
        for i, v in enumerate((0, -1, 4, 5))
    ]
    return rows
