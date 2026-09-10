"""Exact conservative tail and canonical-coefficient intervals."""

from fractions import Fraction
from functools import cache

import sympy as sp
from p8_vacuum_fermion_local_matching import calibration as local
from p8_vacuum_gauge_yukawa_screen import calibration as candidate


def exact(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require a finite exact rational")
    return sp.Rational(value)


def tail_bound(mass, Y):
    m, Y = exact(mass), exact(Y)
    if m < 36 or Y <= 0:
        raise ValueError("Require mF>=36 and Y>0")
    return 300000000 * Y * Y / (144 * m**4)


def canonical_band(tree_b2, tail, slope_lower, slope_upper):
    T, E, lo, hi = map(exact, (tree_b2, tail, slope_lower, slope_upper))
    if T <= 0 or E < 0 or not 0 <= lo <= hi < 1:
        raise ValueError("Require T>0, E>=0 and 0<=slope_lower<=slope_upper<1")
    # A negative numerator reverses denominator monotonicity.
    lower = (T - E) / ((1 - lo) ** 2 if T >= E else (1 - hi) ** 2)
    upper = (T + E) / (1 - hi) ** 2
    return {
        "lower": lower,
        "upper": upper,
        "positive_selected_coefficient": bool(lower > 0),
        "strict_selected_increase_over_tree": bool(lower > T),
        "scope": "Only the selected tree plus fermion functional; not the full quantum candidate.",
    }


@cache
def data():
    p = local.data()
    old = candidate.data()
    m = p["fermion_mass"]
    Yhi = p["rational_Yukawa_squared_upper"]
    Y = p["actual_Yukawa_squared"]
    M = old["old_heavy_mass_squared"]
    g = old["old_cubic_squared"]
    lam = g / (2 * (M - 2) ** 3)
    T = 4 * lam
    E = tail_bound(m, Yhi)
    lo, hi = p["on_shell_slope_bounds"]
    band = canonical_band(T, E, lo, hi)
    formal_lower = 8 * lam * lo - E
    formal_upper = 8 * lam * hi + E
    checks = {
        "same_actual_tree_lambda": lam - sp.Rational(1, 10**600),
        "same_actual_four_scalar_tree_coefficient": T - 4 * lam,
        "same_actual_fermion_mass": m - 10**200,
        "rational_box_tail_bound": E - 300000000 * Yhi**2 / (144 * m**4),
        "formal_field_factor_coefficient": 8 * lam * hi / T - 2 * hi,
        "formal_lower_dictionary": formal_lower + E - 8 * lam * lo,
        "formal_upper_dictionary": formal_upper - E - 8 * lam * hi,
        "exact_normalization_lower_dictionary": band["lower"] * (1 - lo) ** 2 - (T - E),
        "exact_normalization_upper_dictionary": band["upper"] * (1 - hi) ** 2 - (T + E),
    }
    bounds = {
        "actual_Y_below_rational_upper": bool(0 < Y < Yhi),
        "actual_mass_in_complex_tail_domain": bool(m >= 36),
        "box_tail_absolute_below_one_e_minus_1204": bool(
            0 < E < sp.Rational(1, 10**1204)
        ),
        "box_tail_relative_below_one_e_minus_604": bool(
            0 < E / T < sp.Rational(1, 10**604)
        ),
        "positive_small_pole_slopes": bool(0 < lo < hi < sp.Rational(1, 10**205)),
        "formal_fermion_correction_strictly_positive": bool(
            0 < formal_lower < formal_upper
        ),
        "formal_relative_correction_below_one_e_minus_204": bool(
            formal_upper / T < sp.Rational(1, 10**204)
        ),
        "selected_exact_normalized_band_strictly_above_tree": bool(
            T < band["lower"] < band["upper"]
        ),
        "selected_relative_increase_below_one_e_minus_204": bool(
            band["upper"] < T * (1 + sp.Rational(1, 10**204))
        ),
        "selected_coefficient_positive": band["positive_selected_coefficient"],
    }
    return {
        "fermion_mass": m,
        "Yukawa_squared_rational_upper": Yhi,
        "same_reference_tree_lambda": lam,
        "same_reference_tree_second_coefficient": T,
        "rational_box_second_coefficient_absolute_upper": E,
        "rational_box_second_coefficient_relative_upper": E / T,
        "same_fermion_pole_slope_interval": [lo, hi],
        "formal_one_loop_fermion_increment_interval": [formal_lower, formal_upper],
        "formal_one_loop_fermion_relative_increment_upper": formal_upper / T,
        "selected_exact_normalized_coefficient_interval": band,
        "selected_relative_increase_upper": band["upper"] / T - 1,
        "bounds": bounds,
        "scope": "The exact constant-normalized selected tree plus one-loop fermion functional has a second coefficient strictly above its reference tree value, by less than a relative 10^-204. The formal one-loop increment is separately bounded and positive. This is not the full new-model amplitude: old scalar loops, their scheme conversion, later mixed/gauge loops and the truncation error required for a V verdict are not included.",
        "checks": checks,
    }
