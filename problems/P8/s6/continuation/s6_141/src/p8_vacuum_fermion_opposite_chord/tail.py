"""Joint soft projection and exact two-loop radial majorant."""

from functools import cache

import sympy as sp
from p8_vacuum_fermion_self_energy_chord.tail import rational


def enclosure(m, Yhi, ahi, Qlo):
    m, Yhi, ahi, Qlo = map(rational, (m, Yhi, ahi, Qlo))
    if m < 720 or Yhi < 0 or ahi < 0 or not 0 < Qlo <= 144:
        raise ValueError("Need m>=720, nonnegative Y/a and a conservative 0<Qlo<=144")
    pref = sp.Integer(10) ** 14 * 6 / (Qlo**2 * m**4)
    scalar = pref * Yhi**3
    gauge = pref * sp.Rational(16, 3) * ahi * Yhi**2
    return {
        "fermion_mass": m,
        "minimum_joint_Cauchy_radius": m / 360,
        "scalar_opposite_chord_b2_absolute_upper": scalar,
        "gauge_opposite_chord_b2_absolute_upper": gauge,
        "combined_opposite_chord_b2_absolute_upper": scalar + gauge,
        "newly_bounded_words_per_sector": 12,
        "previous_self_energy_words_per_sector": 24,
        "cumulatively_bounded_words_per_sector": 36,
        "remaining_vertex_words_per_sector": 24,
        "scheme": "No divergent proper subgraph. The overall MS quartic contact, including dimensional finite contact terms, has zero b2.",
        "scope": "Only the complete twelve-word opposite-chord subsets; the 24 vertex words in each sector remain unbounded.",
    }


@cache
def data():
    R, x, y, m, t, u, n = sp.symbols("R x y m t u n", positive=True)
    primitive = -2 / sp.sqrt(m * m + x)
    compact = sp.integrate(u * (1 - u), (u, 0, 1))
    radial = 4 * compact
    pref = 12 * 512 * 360**4 * radial
    primitive_n = -2 * (1 + t) ** (-n / 2) / n + 2 * (1 + t) ** (-n / 2 - 1) / (n + 2)
    checks = {
        "geometric_tail_degree_four": sp.factor(
            1 / (R**4 * (1 - 1 / R)) - 1 / (R**3 * (R - 1))
        ),
        "geometric_bound_for_radius_at_least_two": sp.factor(
            2 / R**4 - 1 / (R**3 * (R - 1)) - (R - 2) / (R**4 * (R - 1))
        ),
        "six_fermion_trace_bound": 4 * 2**6 - 256,
        "soft_tail_geometric_factor": 2 * 256 - 512,
        "ordered_angular_radial_measure": sp.factor(x * y / x - y),
        "larger_radius_antiderivative": sp.diff(primitive, x)
        - (m * m + x) ** (-sp.Rational(3, 2)),
        "larger_radius_integral": sp.limit(primitive, x, sp.oo)
        - primitive.subs(x, y)
        - 2 / sp.sqrt(m * m + y),
        "two_radial_regions_after_inner_integral": 2 * 2 - 4,
        "fourth_order_radial_mass_scaling": sp.factor(
            (m * m * t) * m * m / (m * m + m * m * t) ** 4 - t / (m**4 * (1 + t) ** 4)
        ),
        "radial_compactification": sp.factor(
            ((1 - u) / u) / (1 + (1 - u) / u) ** 4 / u**2 - u * (1 - u)
        ),
        "radial_beta_moment": compact - sp.Rational(1, 6),
        "complete_two_loop_radial_moment": radial - sp.Rational(2, 3),
        "positive_soft_degree_antiderivative": sp.simplify(
            sp.diff(primitive_n, t) - t / (1 + t) ** (2 + n / 2)
        ),
        "general_positive_degree_radial_constant": sp.factor(
            -4 * primitive_n.subs(t, 0) - 16 / (n * (n + 2))
        ),
        "degree_four_general_radial_specialization": sp.Rational(16, 4 * 6) - radial,
        "exact_prefactor": pref - 68797071360000,
        "gauge_index_and_open_Casimir_factor": 4 * sp.Rational(4, 3)
        - sp.Rational(16, 3),
        "zero_Yukawa_zero_selected_amplitude": enclosure(720, 0, 1, 144)[
            "combined_opposite_chord_b2_absolute_upper"
        ],
    }
    return {
        "exact_radial_moment": radial,
        "general_positive_soft_degree_moment": "16/[n(n+2)m^n] for n>0",
        "exact_all_word_prefactor": pref,
        "rounded_prefactor": sp.Integer(10) ** 14,
        "rounded_prefactor_is_conservative": bool(pref < 10**14),
        "soft_projection": "Project degrees >=4 before both loop integrals; the full S4/Lorentz subset has degree zero contact, no odd scalar and degree two constant on the equal-mass shell.",
        "forward_coefficient": "The integrated projected amplitude is holomorphic on |s-2|<=1. Cauchy's formula bounds b2 by the same uniform majorant.",
        "scope": "No internal momentum cutoff, imaginary chord shift or assumed vertex-chord subgraph bound.",
        "checks": checks,
    }
