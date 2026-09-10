"""Nonlocal on-shell remainder, explicitly not a finite MS slope certificate."""

from functools import cache

import sympy as sp
from p8_vacuum_fermion_self_energy_chord.tail import rational


def enclosure(m, Y, a, Q):
    m, Y, a, Q = map(rational, (m, Y, a, Q))
    if m < 720 or min(Y, a) < 0 or not 0 < Q <= 144:
        raise ValueError("Need m>=720, nonnegative Y,a and 0<Q<=144")
    pref = 6 * Y / (Q**2 * m**2)
    c = Y + sp.Rational(16, 3) * a
    self_bound = 10**13 * pref * (3 * Y + 16 * a)
    vertex_bound = 3 * 10**13 * pref * c
    total = self_bound + vertex_bound
    return {
        "fermion_mass": m,
        "N": 6,
        "Y_upper": Y,
        "a_upper": a,
        "Q_lower": Q,
        "minimum_Cauchy_radius": m / 360,
        "self_energy_soft_tail_upper": self_bound,
        "vertex_soft_tail_upper": vertex_bound,
        "both_quadratic_primitives_soft_tail_upper": total,
        "outer_s_disc_center": 1,
        "outer_s_disc_radius": 2,
        "inner_on_shell_disc_radius": sp.Rational(1, 2),
        "nonlocal_on_shell_divided_remainder_upper": total / 3,
        "scope": "Only the on-shell-affine-subtracted nonlocal primitive remainder. Finite MS mass and slope, other forests and canonical two-loop pole remain uncomputed.",
    }


@cache
def data():
    s, z, alpha, beta = sp.symbols("s z alpha beta")
    polynomial = alpha + beta * s
    affine = lambda f: sp.expand(f - f.subs(s, 1) - (s - 1) * sp.diff(f, s).subs(s, 1))
    R, r = sp.symbols("outer_radius inner_radius", positive=True)
    return {
        "projection": "The complete two-point function is even and O(4) invariant. Soft degrees zero and two are a mass plus kinetic polynomial; odd degrees vanish. Subtract these before integrations.",
        "on_shell_affine_map": "f_R(s)=f(s)-f(1)-(s-1)f'(1), independent of an arbitrary mass/kinetic polynomial.",
        "second_Cauchy_bound": "If the soft tail is bounded by E on |s-1|<=2, then |f_R(s)|/(s-1)^2<=E/3 on |s-1|<=1/2.",
        "local_reference_boundary": "This identity does not evaluate f_MS(1) or f_MS'(1). It cannot by itself fix the canonical two-loop field factor or certify the complete pole.",
        "checks": {
            "mass_and_kinetic_polynomial_removed": affine(polynomial),
            "arbitrary_local_reference_independence": affine(s**3 + polynomial)
            - affine(s**3),
            "degree_four_soft_becomes_quadratic_invariant": (z**2) ** 2 - z**4,
            "third_invariant_affine_subtraction": sp.factor(
                affine(s**3) - (s - 1) ** 2 * (s + 2)
            ),
            "geometric_Cauchy_remainder_sum": sp.factor(
                1 / (R**2 * (1 - r / R)) - 1 / (R * (R - r))
            ),
            "actual_two_disc_constant": 1
            / (sp.Integer(2) ** 2 * (1 - sp.Rational(1, 2) / 2))
            - sp.Rational(1, 3),
            "no_Yukawa_no_quadratic_increment": enclosure(720, 0, 1, 144)[
                "nonlocal_on_shell_divided_remainder_upper"
            ],
            "combined_coupling_dictionary": (3 * sp.Symbol("Y") + 16 * sp.Symbol("a"))
            - 3 * (sp.Symbol("Y") + sp.Rational(16, 3) * sp.Symbol("a")),
        },
    }
