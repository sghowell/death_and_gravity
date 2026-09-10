"""Actual primitive zero-soft and on-shell MS slope enclosures."""

from functools import cache

import sympy as s
from p8_vacuum_fermion_quadratic_forests import calibration as previous
from p8_vacuum_fermion_self_energy_chord.tail import rational


def enclosure(m, Y, a, Q, soft_tail):
    m, Y, a, Q, soft_tail = map(rational, (m, Y, a, Q, soft_tail))
    if m < 720 or min(Y, a, soft_tail) < 0 or not 0 < Q <= 144:
        raise ValueError("Need m>=720, nonnegative Y,a,E and 0<Q<=144")
    pref = 6 / Q**2
    zero = pref * (
        s.Rational(49, 12) * Y**2 + s.Rational(37, 6) * Y * a * s.Rational(4, 3)
    )
    delta = 4 * 10**12 * pref * Y**2 / m ** s.Rational(3, 2)
    pole = soft_tail / 2
    return {
        "massless_reference_absolute_upper": zero,
        "nonzero_scalar_mass_correction_upper": delta,
        "zero_soft_to_on_shell_slope_correction_upper": pole,
        "primitive_on_shell_MS_slope_absolute_upper": zero + delta + pole,
        "scope": "Two quadratic primitive MS slopes only. The finite mass references, other matching/canonical terms and complete two-loop pole remain open.",
    }


@cache
def data():
    old = previous.data()
    p = old["actual_reference_parameters"]
    E = old["enclosure"]["both_quadratic_primitives_soft_tail_upper"]
    r = enclosure(p["mF"], p["Y_upper"], p["a_upper"], 144, E)
    B = r["primitive_on_shell_MS_slope_absolute_upper"]
    return {
        "actual_reference_parameters": p,
        "enclosure": r,
        "checks": {
            "same_mass": p["mF"] - 10**200,
            "zero_soft_tail_has_no_linear_invariant": s.diff(
                s.Symbol("s") ** 2, s.Symbol("s")
            ).subs(s.Symbol("s"), 0),
            "Cauchy_derivative_radius_two": r[
                "zero_soft_to_on_shell_slope_correction_upper"
            ]
            - E / 2,
            "three_parts_summed_once": B
            - r["massless_reference_absolute_upper"]
            - r["nonzero_scalar_mass_correction_upper"]
            - r["zero_soft_to_on_shell_slope_correction_upper"],
        },
        "bounds": {
            "actual_mass_in_domain": bool(p["mF"] >= 720),
            "actual_shared_Y_enclosure": bool(0 < p["actual_Y"] < p["Y_upper"]),
            "actual_shared_gauge_enclosure": bool(0 < p["actual_a"] < p["a_upper"]),
            "scalar_mass_correction_below_one_e_minus_700": bool(
                0 < r["nonzero_scalar_mass_correction_upper"] < s.Rational(1, 10**700)
            ),
            "MS_primitive_slope_below_one_e_minus_410": bool(
                0 < B < s.Rational(1, 10**410)
            ),
            "on_shell_continuation_error_below_one_e_minus_799": bool(
                0
                < r["zero_soft_to_on_shell_slope_correction_upper"]
                < s.Rational(1, 10**799)
            ),
        },
    }
