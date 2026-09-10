"""Actual local references and mass-one continuation of the single older row."""

from functools import cache

import sympy as s
from p8_polynomial_vacuum import model
from p8_vacuum_fermion_ms_mass import calibration as previous

from . import bubble, tadpole


@cache
def data():
    p = previous.data()["actual_reference_parameters"]
    v = model.data()["actual_parameters"]
    m, Y = p["mF"], p["Y_upper"]
    g, M, L = (
        v["cubic_coupling_squared"],
        v["heavy_mass_squared"],
        v["bare_polynomial_quartic"],
    )
    G = s.Rational(1, 8192)
    t = tadpole.enclosure(m, Y, 144)
    b = bubble.enclosure(m, Y, M, 144)
    local = L * t["complete_tadpole_absolute_upper"] / 2
    mixed = g * b["complete_finite_bubble_absolute_upper"]
    zero = local + mixed
    slope = g * b["uniform_derivative_per_unit_g_upper"]
    source = G * t["complete_tadpole_absolute_upper"] / 2
    return {
        "actual_reference_parameters": {**p, "g": g, "M": M, "L": L, "G": G},
        "finite_inserted_tadpole_enclosure": t,
        "finite_heavy_bubble_enclosure": b,
        "local_tadpole_mass_absolute_upper": local,
        "heavy_bubble_mass_absolute_upper": mixed,
        "complete_row_zero_mass_reference_absolute_upper": zero,
        "uniform_row_MS_slope_absolute_upper": slope,
        "complete_row_on_shell_mass_reference_absolute_upper": zero + slope,
        "assigned_H_source_absolute_upper": source,
        "assigned_H_shift_absolute_upper": source / M,
        "checks": {
            "same_actual_fermion_mass": m - 10**200,
            "same_actual_cubic_square": G * G - g,
            "both_local_mass_pieces_retained": zero - local - mixed,
            "same_ratio_dyadic": b["dyadic_T_over_M"] - 676,
            "fixed_H_source_factor": 2 * source
            - G * t["complete_tadpole_absolute_upper"],
        },
        "bounds": {
            "actual_scalar_heavy_spectral_hierarchy": bool(1 < M < 4 * m * m / 16),
            "actual_positive_Y_enclosure": bool(0 < p["actual_Y"] < Y),
            "actual_leading_couplings_positive": bool(min(L, g, G) > 0),
            "finite_tadpole_remainder_below_one_e_minus_606": bool(
                0
                < t["finite_tadpole_remainder_absolute_upper"]
                < s.Rational(1, 10**606)
            ),
            "finite_bubble_mass_below_one_e_minus_214": bool(
                0 < mixed < s.Rational(1, 10**214)
            ),
            "complete_row_mass_below_one_e_minus_10": bool(
                0 < zero + slope < s.Rational(1, 10**10)
            ),
            "assigned_H_shift_below_one_e_minus_6": bool(
                0 < source / M < s.Rational(1, 10**6)
            ),
            "row_slope_below_one_e_minus_616": bool(0 < slope < s.Rational(1, 10**616)),
        },
        "scope": "Local references of the older paired insertion row only. The two vacuum rows, other matching/counterterm and canonical terms, full pole/error and original P8 remain open.",
    }
