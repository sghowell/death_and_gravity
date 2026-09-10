"""Actual complete primitive MS mass references, with on-shell continuation."""

from functools import cache

import sympy as s
from p8_vacuum_fermion_ms_slopes import calibration as slope
from p8_vacuum_fermion_quadratic_forests import calibration as nonlocal_part

from . import correction


@cache
def data():
    previous = slope.data()
    p = previous["actual_reference_parameters"]
    m, Y, a = p["mF"], p["Y_upper"], p["a_upper"]
    mass_diff = correction.enclosure(m, Y, 144)
    pref = s.Rational(6, 144**2)
    scalar0 = pref * 56 * Y**2 * m * m
    gauge0 = pref * 40 * Y * a * s.Rational(4, 3) * m * m
    zero = (
        scalar0 + gauge0 + mass_diff["complete_scalar_mass_difference_absolute_upper"]
    )
    z = previous["enclosure"]
    zero_slope = (
        z["massless_reference_absolute_upper"]
        + z["nonzero_scalar_mass_correction_upper"]
    )
    E = nonlocal_part.data()["enclosure"]["both_quadratic_primitives_soft_tail_upper"]
    on_shell = zero + zero_slope + E
    return {
        "actual_reference_parameters": p,
        "finite_scalar_mass_difference": mass_diff,
        "massless_scalar_reference_absolute_upper": scalar0,
        "massless_gauge_reference_absolute_upper": gauge0,
        "actual_zero_momentum_MS_mass_reference_absolute_upper": zero,
        "actual_zero_soft_MS_slope_absolute_upper": zero_slope,
        "previous_soft_tail_upper": E,
        "actual_on_shell_MS_mass_reference_absolute_upper": on_shell,
        "checks": {
            "same_mass": m - 10**200,
            "actual_zero_reference_all_pieces": zero
            - scalar0
            - gauge0
            - mass_diff["complete_scalar_mass_difference_absolute_upper"],
            "on_shell_mass_reference_decomposition": on_shell - zero - zero_slope - E,
            "same_threshold_dyadic": mass_diff["dyadic_threshold"] - 1331,
            "same_mass_squared_dyadic": mass_diff["dyadic_m_squared"] - 1329,
        },
        "bounds": {
            "actual_mass_in_domain": bool(m >= 720),
            "actual_positive_Y_enclosure": bool(0 < p["actual_Y"] < Y),
            "actual_positive_gauge_enclosure": bool(0 < p["actual_a"] < a),
            "scalar_mass_difference_below_one_e_minus_408": bool(
                0
                < mass_diff["complete_scalar_mass_difference_absolute_upper"]
                < s.Rational(1, 10**408)
            ),
            "finite_ratio_remainder_below_one_e_minus_808": bool(
                0
                < mass_diff["finite_mass_ratio_remainder_absolute_upper"]
                < s.Rational(1, 10**808)
            ),
            "zero_momentum_reference_below_one_e_minus_10": bool(
                0 < zero < s.Rational(1, 10**10)
            ),
            "on_shell_reference_below_one_e_minus_10": bool(
                0 < on_shell < s.Rational(1, 10**10)
            ),
        },
        "scope": "The two quadratic primitive finite MS mass references and prior slopes/nonlocal parts are controlled. This does not combine other matching forests, the older W1F0 sector, vacuum rows or the full canonical two-loop pole.",
    }
