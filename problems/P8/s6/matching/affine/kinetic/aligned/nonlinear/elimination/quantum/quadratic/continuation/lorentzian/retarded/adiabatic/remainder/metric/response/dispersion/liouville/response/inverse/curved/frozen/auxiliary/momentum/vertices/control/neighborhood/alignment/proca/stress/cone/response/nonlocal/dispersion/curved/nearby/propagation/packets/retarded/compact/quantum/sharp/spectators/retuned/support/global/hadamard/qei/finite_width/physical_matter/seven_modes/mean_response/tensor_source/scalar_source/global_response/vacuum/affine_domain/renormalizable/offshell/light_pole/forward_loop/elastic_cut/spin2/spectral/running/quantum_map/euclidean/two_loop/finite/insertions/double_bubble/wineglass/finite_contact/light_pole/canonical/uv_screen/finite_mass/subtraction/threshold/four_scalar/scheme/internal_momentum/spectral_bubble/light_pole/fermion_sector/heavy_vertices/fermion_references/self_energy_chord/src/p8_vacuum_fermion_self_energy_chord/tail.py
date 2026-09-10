"""Soft Taylor tail and exact radial integration of the self-energy chords."""

from functools import cache

import sympy as sp


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, sp.Integer, sp.Rational)):
        raise TypeError("Require an exact finite rational")
    return sp.Rational(value)


def enclosure(m, Yhi, ahi, Qlo):
    m, Yhi, ahi, Qlo = map(rational, (m, Yhi, ahi, Qlo))
    if m < 720 or Yhi < 0 or ahi < 0 or not 0 < Qlo <= 144:
        raise ValueError("Need m>=720, nonnegative Y/a and a conservative 0<Qlo<=144")
    pref = sp.Integer(2) * 10**14 * 6 / (Qlo**2 * m**4)
    scalar = pref * 3 * Yhi**3
    gauge = pref * 16 * ahi * Yhi**2
    return {
        "fermion_mass": m,
        "minimum_soft_scaling_Cauchy_radius": m / 360,
        "scalar_self_energy_chord_b2_absolute_upper": scalar,
        "gauge_self_energy_chord_b2_absolute_upper": gauge,
        "combined_self_energy_chord_b2_absolute_upper": scalar + gauge,
        "selected_words_per_sector": 24,
        "unbounded_words_per_sector": 36,
        "scheme": "Complete proper MS self-energy and its counterterms, using leading shared reference couplings; the remaining overall quartic contact has zero b2.",
        "scope": "Uniform projected-amplitude and forward b2 bounds for the 24 self-energy chords in each sector only.",
    }


@cache
def data():
    R = sp.Symbol("Cauchy_radius", positive=True)
    u, t, m = sp.symbols("u t m", positive=True)
    I0 = sp.integrate(u * (1 - u), (u, 0, 1))
    Ilog = sp.integrate(-u * (1 - u) * sp.log(u), (u, 0, 1))
    radial = Ilog + 6 * I0
    pref = 24 * 256 * 360**4 * radial
    checks = {
        "geometric_tail_from_degree_four": sp.factor(
            1 / (R**4 * (1 - 1 / R)) - 1 / (R**3 * (R - 1))
        ),
        "geometric_tail_below_twice_fourth_power": sp.factor(
            2 / R**4 - 1 / (R**3 * (R - 1)) - (R - 2) / (R**4 * (R - 1))
        ),
        "five_shifted_fermion_propagators_and_Dirac_trace": 4 * 2**5 - 128,
        "geometric_tail_factor_two": 2 * 128 - 256,
        "all_twenty_four_marked_box_words": 6 * 4 - 24,
        "radial_mass_scaling": sp.factor(
            (m * m * t) * m * m / (m * m + m * m * t) ** 4 - t / (m**4 * (1 + t) ** 4)
        ),
        "radial_compactification_density": sp.factor(
            ((1 - u) / u) / (1 + (1 - u) / u) ** 4 / u**2 - u * (1 - u)
        ),
        "radial_without_log": I0 - sp.Rational(1, 6),
        "radial_logarithmic_moment": Ilog - sp.Rational(5, 36),
        "complete_radial_moment": radial - sp.Rational(41, 36),
        "exact_combinatoric_radial_prefactor": pref
        - sp.Integer(24) * 256 * 360**4 * sp.Rational(41, 36),
        "active_closed_gauge_color_factor": 6 * sp.Rational(4, 3) - 8,
        "open_Casimir_in_gauge_majorant": 12 * sp.Rational(4, 3) - 16,
        "no_Yukawa_no_scalar_or_gauge_box": enclosure(720, 0, 1, 144)[
            "combined_self_energy_chord_b2_absolute_upper"
        ],
    }
    return {
        "exact_radial_moment": radial,
        "exact_all_word_prefactor": pref,
        "rounded_prefactor": sp.Integer(2) * 10**14,
        "rounded_prefactor_is_conservative": bool(pref < 2 * 10**14),
        "soft_projection": "Remove degrees zero through three before outer integration. Complete-subset Lorentz/S4 symmetry makes degree zero a contact, odd degrees vanish, and degree two proportional to sum external p_i^2, hence constant on the mass-one shell.",
        "Cauchy_steps": "First bound the soft-scaling tail pointwise with R(q)=sqrt(mF^2+q^2)/360; then integrate its uniform majorant. The resulting holomorphic unit-forward-disc bound gives the same bound on b2.",
        "scope": "No q/m expansion, finite momentum cutoff, unpaired UV integral, or transfer of this Cauchy radius to the remaining vertex-chord graphs.",
        "checks": checks,
    }
