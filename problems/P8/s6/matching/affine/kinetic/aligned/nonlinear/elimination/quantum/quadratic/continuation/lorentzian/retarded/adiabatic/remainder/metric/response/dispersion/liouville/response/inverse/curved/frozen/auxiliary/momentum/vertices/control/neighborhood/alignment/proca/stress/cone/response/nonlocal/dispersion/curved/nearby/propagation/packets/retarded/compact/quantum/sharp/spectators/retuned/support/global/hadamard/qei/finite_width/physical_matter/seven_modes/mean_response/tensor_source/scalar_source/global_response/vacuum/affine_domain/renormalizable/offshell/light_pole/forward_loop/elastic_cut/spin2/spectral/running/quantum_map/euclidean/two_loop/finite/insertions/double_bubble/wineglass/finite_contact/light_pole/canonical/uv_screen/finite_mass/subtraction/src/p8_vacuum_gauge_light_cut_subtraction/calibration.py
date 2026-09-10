"""Rational coefficient enclosure for the unchanged candidate's finite cut."""

from fractions import Fraction
from functools import cache

import sympy as sp
from p8_vacuum_finite_mass_gauge_cut import calibration as finite_mass
from p8_vacuum_gauge_yukawa_screen import calibration as leading

from . import logarithm

DEFAULT_MASS = sp.Integer(10) ** 200


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require a finite exact rational fermion mass")
    return sp.Rational(value)


def point(mass=DEFAULT_MASS):
    m = rational(mass)
    if m < 36:
        raise ValueError("Require mF>=36 for the complex routing domain")
    eps = 4000000 / m
    eta = sp.Rational(5, 2) * eps + eps * eps / 8
    error = 288 * eta
    log = logarithm.enclosure()
    low = 2 * log["lower"] - sp.Rational(21, 4) - error
    high = 2 * log["upper"] - sp.Rational(21, 4) + error
    return {
        "fermion_mass": m,
        "complex_Neumann_ratio_upper": 18 / m,
        "box_epsilon": eps,
        "normalized_rho_error_eta": eta,
        "finite_mass_cut_b2_error_over_K_upper": error,
        "cut_b2_over_K_lower": low,
        "cut_b2_over_K_upper": high,
        "strictly_negative_subtraction_coefficient": bool(high < 0),
        "strict_minus_four_to_minus_three_K_band": bool(-4 < low < high < -3),
        "scope": "Only the coefficient of the prescribed first-two-gauge finite-cut subtraction; not the candidate's full b2 or a positivity verdict.",
    }


@cache
def data():
    old = finite_mass.data()
    pre = leading.data()
    m = old["fermion_mass"]
    z = old["box_scale_z"]
    K = 8 * z * z / (9 * sp.pi**2)
    upper = pre["leading_gauge_bubble_log_coefficient_rational_upper"]
    d = point(m)
    checks = {
        "same_fermion_mass": m - 10**200,
        "same_K_from_inherited_operator": sp.simplify(
            K
            - 8
            * pre["leading_threshold_coefficient_at_prospective_boundary"] ** 2
            / sp.pi**2
        ),
        "actual_complex_Neumann_ratio": d["complex_Neumann_ratio_upper"]
        - sp.Rational(18, 10**200),
        "actual_complex_box_epsilon": d["box_epsilon"] - sp.Rational(4, 10**194),
        "actual_rho_error": d["normalized_rho_error_eta"]
        - (sp.Rational(10, 10**194) + sp.Rational(2, 10**388)),
        "actual_coefficient_error_dictionary": d[
            "finite_mass_cut_b2_error_over_K_upper"
        ]
        - 288 * d["normalized_rho_error_eta"],
        "coefficient_interval_width": d["cut_b2_over_K_upper"]
        - d["cut_b2_over_K_lower"]
        - 2 * logarithm.enclosure()["tail_upper"]
        - 2 * d["finite_mass_cut_b2_error_over_K_upper"],
    }
    bounds = {
        "same_positive_finite_mass_in_complex_domain": m >= 36,
        "actual_strict_complex_Neumann_disc": 0 < 18 / m < sp.Rational(1, 2),
        "positive_cut_normalization": K > 0,
        "actual_coefficient_error_below_one_e_minus_190_K": 0
        < d["finite_mass_cut_b2_error_over_K_upper"]
        < sp.Rational(1, 10**190),
        "strict_negative_subtraction_coefficient": d[
            "strictly_negative_subtraction_coefficient"
        ],
        "strict_normalized_band_minus_four_to_minus_three": d[
            "strict_minus_four_to_minus_three_K_band"
        ],
        "rational_K_upper_positive_below_one_e_minus_1620": 0
        < upper
        < sp.Rational(1, 10**1620),
        "absolute_cut_coefficient_below_four_e_minus_1620": 4 * upper
        < sp.Rational(4, 10**1620),
        "small_mass_control_not_forced_to_a_sign": not point(36)[
            "strictly_negative_subtraction_coefficient"
        ],
    }
    return {
        "candidate": "GY14-unbroken with unchanged S6.128/S6.129 prospective boundary",
        "K": K,
        "K_rational_upper": upper,
        "absolute_subtraction_b2_rational_upper": 4 * upper,
        "actual": d,
        "small_mass_control": point(36),
        "checks": checks,
        "bounds": bounds,
        "scope": "The finite-mass subtraction coefficient, including both known cuts, is enclosed in (-4K,-3K). This is not a new-model amplitude budget, an adjustable renormalizable counterterm, or a proof of high-energy boundedness or original V/G/B closure.",
    }


def bad_cases():
    invalid = (
        True,
        False,
        1.0,
        sp.Float(36),
        "36",
        None,
        sp.oo,
        -sp.oo,
        sp.zoo,
        sp.I,
        sp.nan,
        sp.Symbol("m"),
    )
    return [(f"inexact_mass_{j}", point, (v,)) for j, v in enumerate(invalid)] + [
        (f"outside_complex_mass_domain_{j}", point, (v,))
        for j, v in enumerate((-1, 0, 35))
    ]
