"""Rational full-family enclosure from absolutely convergent spectral integrals."""

from functools import cache

import sympy as sp


def rational(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, sp.Integer, sp.Rational)):
        raise TypeError(name + " must be an exact rational")
    return sp.Rational(value)


def dyadic_log_upper(value):
    value = rational(value, "logarithm argument")
    if value < 1:
        raise ValueError("logarithm argument must be at least one")
    n = max(0, int(value.p).bit_length() - int(value.q).bit_length())
    while sp.Integer(2) ** n < value:
        n += 1
    while n and sp.Integer(2) ** (n - 1) >= value:
        n -= 1
    return sp.Integer(n)


def enclosure(m, Y, L, g, M, Qlo):
    m, Y, L, g, M, Qlo = (
        rational(v, name)
        for v, name in zip((m, Y, L, g, M, Qlo), ("m", "Y", "L", "g", "M", "Qlo"))
    )
    if m < 2 or Y < 0 or L < 0 or g < 0 or M < 24 or Qlo <= 0:
        raise ValueError("need m>=2, nonnegative Y,L,g, M>=24 and Qlo>0")
    T = 4 * m * m
    C = 12 * Y / Qlo
    V = L + g / (M - 3)
    n = dyadic_log_upper(4 * T)
    ell = n + 1
    bubble = 36 * C * V**2 / (Qlo * T)
    triangle = 96 * C * g * V * ell / (Qlo * T)
    box = 192 * C * g * g / (Qlo * T)
    return {
        "spectral_threshold": T,
        "spectral_density_C_upper": C,
        "external_vertex_absolute_upper": V,
        "dyadic_exponent": n,
        "log_4T_plus_one_upper": ell,
        "all_channel_bubble_amplitude_upper": bubble,
        "all_channel_triangle_amplitude_upper": triangle,
        "all_channel_box_amplitude_upper": box,
        "full_family_disc_amplitude_and_b2_upper": bubble + triangle + box,
        "scheme": "paired inner physical on-shell references; outer local full-parent zero-bubble reference",
        "not_included": "other primitive families, other counterterms and common finite MS/canonical conversion",
    }


@cache
def data():
    v, T, C, Q, V, g = sp.symbols("v T C Q V g", positive=True)
    logprimitive = -(sp.log(4 * v) + 1) / v
    bubble_tail = 12 * C / T
    triangle_tail = 8 * C * (sp.log(4 * T) + 1) / T
    box_tail = 16 * C / T
    x = sp.Symbol("x", nonnegative=True)
    return {
        "spectral_bubble_absolute_upper": bubble_tail,
        "spectral_triangle_absolute_upper": triangle_tail,
        "spectral_box_absolute_upper": box_tail,
        "full_family_upper": C
        / (Q * T)
        * (36 * V * V + 96 * g * V * (sp.log(4 * T) + 1) + 192 * g * g),
        "checks": {
            "log_weighted_tail_primitive": sp.simplify(
                sp.diff(logprimitive, v) - sp.log(4 * v) / v**2
            ),
            "log_weighted_tail_anchor": sp.simplify(
                -logprimitive.subs(v, T) - (sp.log(4 * T) + 1) / T
            ),
            "inverse_mass_tail_primitive": sp.diff(-1 / v, v) - 1 / v**2,
            "inverse_mass_tail_anchor": -(-1 / v).subs(v, T) - 1 / T,
            "subtracted_bubble_channel_bound": bubble_tail
            - 3 * sp.Rational(1, 2) * 8 * C / T,
            "triangle_spectral_tail_constant": triangle_tail
            - 4 * C * 2 * (sp.log(4 * T) + 1) / T,
            "box_spectral_tail_constant": box_tail - 4 * C * 4 / T,
            "full_bubble_three_channels": 3 * V**2 * bubble_tail / Q
            - 36 * C * V**2 / (Q * T),
            "four_triangles_three_channels": 3 * 4 * g * V * triangle_tail / Q
            - 96 * C * g * V * (sp.log(4 * T) + 1) / (Q * T),
            "four_boxes_three_channels": 3 * 4 * g**2 * box_tail / Q
            - 192 * C * g**2 / (Q * T),
            "log_two_integral_is_below_one": 1 - 1 / (1 + x) - x / (1 + x),
            "dyadic_at_one": dyadic_log_upper(1),
            "dyadic_at_power_two": dyadic_log_upper(16) - 4,
            "dyadic_just_above_power_two": dyadic_log_upper(sp.Rational(33, 2)) - 5,
            "no_fermion_coupling_no_family": enclosure(2, 0, 1, 1, 24, 144)[
                "full_family_disc_amplitude_and_b2_upper"
            ],
        },
        "scope": "Uniform absolute amplitude bound on the unit disc implies the same b2 bound by Cauchy. No sign for the complete family is inferred from the positive insertion density.",
    }
