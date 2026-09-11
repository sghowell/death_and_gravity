"""Actual SAT8 mass floor and its conditional clock ratios."""

from functools import cache

import sympy as s
from p8_exceptional_vacuum import analytic
from p8_polynomial_vacuum import model as scalar_model
from p8_vacuum_fermion_local_matching import calibration as fermion
from p8_vacuum_gauge_yukawa_screen import flow

from . import clock, profile


@cache
def data():
    p = fermion.data()
    m, Y = p["fermion_mass"], p["rational_Yukawa_squared_upper"]
    R = s.Integer(10) ** 300
    eta = s.Rational(1, 100)
    gap = profile.gap(m, Y, R, eta)
    # A slightly sharper rational cap proves the strict first-power gap:
    # Y<9e-206 implies y<3e-103.
    ycap = s.Rational(3, 10**103)
    ratios = clock.ratios(m, ycap, R, s.sqrt(analytic.KAPPA), 0, eta)
    scalar = scalar_model.data()["actual_parameters"]
    Ylower = s.Rational(19, 90) * scalar["bare_polynomial_quartic"]
    ray = flow.data()["quartic_to_gauge_squared"]
    ylower = s.Rational(2, 10**103)
    lam, gamma = scalar["lambda"], scalar["gamma"]
    map_slope = 1 + 2 * lam * analytic.KAPPA - 2 * gamma * analytic.KAPPA
    naive_ratio_lower = ylower * map_slope * s.sqrt(analytic.KAPPA) / m**2
    return {
        "candidate_name": "GY14-SAT8",
        "fixed_profile_radius": R,
        "same_fermion_mass": m,
        "actual_Yukawa_squared_rational_upper": Y,
        "first_power_Yukawa_rational_upper": ycap,
        "strict_pointwise_paired_mass_enclosure": gap,
        "diagnostic_linear_argument_mass_derivative_ratios": ratios,
        "flat_linear_clock_literal_map_slope_at_zero": map_slope,
        "flat_linear_clock_literal_map_mass_derivative_ratio_lower": naive_ratio_lower,
        "naive_map_screen": "Even on the flat diagnostic Psi=sqrt(kappa)t, the literal cubic field map has slope 1+2lambda kappa-2gamma kappa at zero. The profile has f'(0)=1, so its first mass-derivative ratio there exceeds 1e97. The small direct-argument ratios cannot be transferred to this extrapolated map. This fails a small-derivative expansion screen, not a proof of particle production or exclusion of SAT8, its parent or a ladder row.",
        "bounds": {
            "same_positive_actual_fermion_mass": bool(m == 10**200),
            "squared_Yukawa_cap_proves_first_power_cap": bool(Y < ycap * ycap),
            "actual_global_relative_mass_shift_below_one_percent": bool(
                Y * R * R / m**2 < eta**2
            ),
            "actual_mass_floor_above_ninety_nine_percent": bool(
                gap["both_mass_strict_lower"] == s.Rational(99, 100) * m
            ),
            "diagnostic_first_mass_derivative_ratio_below_one_e_minus_102": bool(
                ratios["first_mass_derivative_over_mass_squared_upper"]
                < s.Rational(1, 10**102)
            ),
            "diagnostic_second_mass_derivative_ratio_below_one_e_minus_201": bool(
                ratios["second_mass_derivative_over_mass_cubed_upper"]
                < s.Rational(1, 10**201)
            ),
            "saturation_radius_below_unprotected_mass_crossing": bool(
                Y * R * R < m * m
            ),
            "positive_ray_rational_lower_cap_proves_Yukawa_lower": bool(
                Ylower > ylower * ylower
            ),
            "same_actual_positive_ray_strictly_below_two": bool(0 < ray < 2),
            "flat_linear_clock_naive_map_mass_smallness_screen_fails": bool(
                naive_ratio_lower > 10**97
            ),
        },
        "checks": {
            "profile_radius_prescription": R - 10**300,
            "same_original_fermion_mass": m - fermion.data()["fermion_mass"],
            "squared_relative_mass_shift": gap["relative_mass_shift_squared_upper"]
            - Y * R * R / m**2,
            "strict_paired_mass_floor": gap["both_mass_strict_lower"] - (1 - eta) * m,
            "strict_paired_mass_ceiling": gap["both_mass_strict_upper"] - (1 + eta) * m,
            "diagnostic_clock_velocity": s.sqrt(analytic.KAPPA) - 10**400,
            "naive_map_slope_on_same_actual_coefficients": map_slope
            - 1
            - 2 * lam * analytic.KAPPA
            + 2 * gamma * analytic.KAPPA,
            "naive_map_ratio_uses_mass_at_zero_not_global_lower": naive_ratio_lower
            * m**2
            - ylower * map_slope * s.sqrt(analytic.KAPPA),
            "lower_Yukawa_squared_from_same_positive_boundary": Ylower
            - s.Rational(19, 45) * scalar["bare_polynomial_quartic"] / 2,
        },
    }
