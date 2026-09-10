"""Exact regional moments and complete paired vertex-chord bounds."""

from functools import cache

import sympy as sp
from p8_vacuum_fermion_self_energy_chord.tail import rational


def enclosure(m, Yhi, ahi, Qlo):
    m, Yhi, ahi, Qlo = map(rational, (m, Yhi, ahi, Qlo))
    if m < 720 or Yhi < 0 or ahi < 0 or not 0 < Qlo <= 144:
        raise ValueError("Need m>=720, nonnegative Y/a and conservative 0<Qlo<=144")
    pref = sp.Integer(5) * 10**14 * 6 / (Qlo**2 * m**4)
    scalar = pref * Yhi**3
    gauge = pref * sp.Rational(16, 3) * ahi * Yhi**2
    return {
        "fermion_mass": m,
        "minimum_Cauchy_radius": m / 360,
        "scalar_vertex_chord_b2_absolute_upper": scalar,
        "gauge_vertex_chord_b2_absolute_upper": gauge,
        "combined_vertex_chord_b2_absolute_upper": scalar + gauge,
        "newly_bounded_words_per_sector": 24,
        "cumulatively_bounded_words_per_sector": 60,
        "remaining_words_in_these_two_primitive_rows": 0,
        "scheme": "The complete proper MS Yukawa vertex is the convergent subtracted kernel plus its finite zero-momentum MS anchor. Overall quartic contacts have zero b2.",
        "scope": "Only two complete quartic primitive rows after adding S6.140/S6.141; not all primitive rows or the full matched two-loop amplitude.",
    }


@cache
def data():
    x, y, m, u, R = sp.symbols("x y m u R", positive=True)
    I0 = sp.integrate(u * (1 - u), (u, 0, 1))
    Ilog = sp.integrate(-u * (1 - u) * sp.log(u), (u, 0, 1))
    raw = 3 * I0
    sub = Ilog + 2 * I0
    high = 32 * I0
    difference = 4 * raw + sub + high
    total = difference + 2 * I0
    base = sp.Integer(24) * 128 * 360**4
    checks = {
        "geometric_soft_tail": sp.factor(
            1 / (R**4 * (1 - 1 / R)) - 1 / (R**3 * (R - 1))
        ),
        "geometric_tail_upper_gap": sp.factor(
            2 / R**4 - 1 / (R**3 * (R - 1)) - (R - 2) / (R**4 * (R - 1))
        ),
        "raw_low_order_y_below_x_inner_integral": sp.integrate(
            1 / (m * m + x) ** 2, (x, y, sp.oo)
        )
        - 1 / (m * m + y),
        "raw_low_first_ordered_region": I0 - sp.Rational(1, 6),
        "raw_low_second_ordered_region": 2 * I0 - sp.Rational(1, 3),
        "raw_low_total_radial_upper": raw - sp.Rational(1, 2),
        "subtraction_low_log_moment": sub - sp.Rational(17, 36),
        "high_radial_inner_integral": sp.integrate(
            y ** (-sp.Rational(3, 2)), (y, 4 * (m * m + x), sp.oo)
        )
        - 1 / sp.sqrt(m * m + x),
        "high_paired_radial_constant": high - sp.Rational(16, 3),
        "difference_all_region_constant": difference - sp.Rational(281, 36),
        "MS_local_anchor_radial_moment": I0 - sp.Rational(1, 6),
        "full_vertex_constant_after_local_anchor_bound": total - sp.Rational(293, 36),
        "exact_all_word_base": base - 51597803520000,
        "gauge_caps_Casimir": 4 * sp.Rational(4, 3) - sp.Rational(16, 3),
        "zero_Yukawa_no_selected_contribution": enclosure(720, 0, 1, 144)[
            "combined_vertex_chord_b2_absolute_upper"
        ],
    }
    return {
        "raw_LOW_radial_upper": raw,
        "subtraction_LOW_radial_upper": sub,
        "paired_HIGH_radial_upper_including_32": high,
        "combined_difference_constant": difference,
        "MS_local_anchor_radial_moment": I0,
        "complete_constant_with_anchor_majorized": total,
        "exact_all_word_base": base,
        "exact_prefactor": base * total,
        "rounded_prefactor": sp.Integer(5) * 10**14,
        "rounded_prefactor_is_conservative": bool(base * total < 5 * 10**14),
        "positive_soft_degree_control": "For every n>0 all regional bounds are integrals of t/(1+t)^(2+n/2), optionally times log(1+t), and converge.",
        "scope": "The full S4/Lorentz paired subset has no b2 below soft degree four; Cauchy projection precedes both integrations.",
        "checks": checks,
    }
