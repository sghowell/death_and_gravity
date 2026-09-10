"""Prospective boundary calibration; no silent matching to the old amplitude."""

from functools import cache

import sympy as sp
from p8_polynomial_vacuum import model as scalar_model
from p8_vacuum_two_loop_canonical import calibration as canonical

from . import flow


@cache
def data():
    old = scalar_model.data()["actual_parameters"]
    L = old["bare_polynomial_quartic"]
    M = old["heavy_mass_squared"]
    g = old["cubic_coupling_squared"]
    rL = flow.data()["quartic_to_gauge_squared"]
    ry = flow.data()["Yukawa_squared_to_gauge_squared"]
    a = L / rL
    Y = ry * a
    mF = sp.Integer(10) ** 200
    C = -a * Y / (48 * sp.pi**2 * mF**2)
    C_upper = L * L / (432 * mF**2)
    weakening = 22 * L * 600 / 144
    log_scale_upper = 600 - 144 / (22 * L)
    cut_coefficient_upper = 8 * C_upper**2 / 9
    old_budget = canonical.point(2)
    checks = {
        "same_prospective_quartic_boundary": sp.simplify(rL * a - L),
        "same_prospective_Yukawa_ray": sp.simplify(Y - ry * a),
        "rational_threshold_majorant_dictionary": C_upper - L * L / (48 * 3**2 * mF**2),
        "rational_running_window_dictionary": weakening
        - 2 * 11 * L * (200 * 3) / (16 * 3**2),
        "rational_log_transmutation_dictionary": log_scale_upper
        - (200 * 3 - 16 * 3**2 / (2 * 11 * L)),
    }
    return {
        "old_scalar_quartic_used_only_as_prospective_boundary": L,
        "old_heavy_mass_squared": M,
        "old_cubic_squared": g,
        "new_common_Dirac_mass": mF,
        "prospective_gauge_squared": a,
        "prospective_Yukawa_squared": Y,
        "leading_threshold_coefficient_at_prospective_boundary": C,
        "leading_threshold_absolute_rational_upper": C_upper,
        "leading_gauge_bubble_log_coefficient_rational_upper": cut_coefficient_upper,
        "one_loop_running_relative_denominator_change_upper_to_light_mass": weakening,
        "one_loop_log_transmutation_scale_upper": log_scale_upper,
        "old_canonical_amplitude_budget_NOT_new_model_amplitude": old_budget,
        "checks": checks,
        "bounds": {
            "prospective_scalar_quartic_positive": L > 0,
            "prospective_scalar_quartic_below_three_e_minus_205": L
            < sp.Rational(3, 10**205),
            "positive_quartic_ray_between_one_and_two": bool(1 < rL < 2),
            "positive_Yukawa_ray_below_one": 0 < ry < 1,
            "new_fermions_above_old_heavy_scalar": mF * mF > M,
            "new_fermions_below_named_Planck_scale": mF < 10**400,
            "leading_operator_rational_upper_below_one_e_minus_810": C_upper
            < sp.Rational(1, 10**810),
            "leading_cut_rational_coefficient_below_one_e_minus_1620": cut_coefficient_upper
            < sp.Rational(1, 10**1620),
            "one_loop_gauge_running_stays_weak_to_light_scale": 0
            < weakening
            < sp.Rational(1, 10**200),
            "formal_transmutation_log_below_minus_one_e_204": log_scale_upper
            < -(10**204),
            "old_scalar_two_loop_budget_positive": old_budget[
                "strictly_positive_canonical_b2_at_this_order"
            ],
        },
        "scope": "L at the new fermion threshold is assigned the old scalar numerical quartic only to screen a possible boundary condition. Finite threshold counterterms and changes to the old amplitude are NOT computed. The transmutation scale is a one-loop running parameter, not a proven Yang-Mills confinement or glueball mass.",
    }
