"""Literal local-parent parameter map from the full subtraction to MS."""

from functools import cache

import sympy as sp


@cache
def data():
    L, g, M, F, G = sp.symbols("L g M F_MS G")
    s, t, u = sp.symbols("s t u")
    channels = (s, t, u)
    C = [-L + g / (M - z) for z in channels]
    A = -L + sum(g / (M - z) for z in channels)
    dL, dg, dM = -3 * L * L * F, -2 * L * g * F, -g * F
    dA = sp.diff(A, L) * dL + sp.diff(A, g) * dg + sp.diff(A, M) * dM
    tree = 2 * g / (M - 2) ** 3
    change = sp.diff(tree, g) * dg + sp.diff(tree, M) * dM
    coeff = -4 * L * g / (M - 2) ** 3 + 6 * g * g / (M - 2) ** 4
    h = sp.Symbol("loop_marker")
    return {
        "matched_star_parameter_shifts": {
            "delta_L": dL,
            "delta_g": dg,
            "delta_M": dM,
            "delta_G": -L * G * F,
        },
        "MS_minus_subtracted_amplitude": "F_MS sum_z[-L+g/(M-z)]^2",
        "MS_minus_subtracted_b2": F * coeff,
        "b2_conversion_coefficient": coeff,
        "relative_to_tree_conversion_coefficient": -2 * L + 3 * g / (M - 2),
        "scope": "Only the order-two outer interaction-reference conversion. Lower-order squared shifts and other canonical-field/parameter insertions remain separately owned.",
        "checks": {
            "literal_full_tree_parameter_variation": sp.factor(
                dA - F * sum(c * c for c in C)
            ),
            "forward_b2_parameter_variation": sp.factor(change - F * coeff),
            "same_frozen_reference_change_coefficient": sp.factor(
                coeff / tree - (-2 * L + 3 * g / (M - 2))
            ),
            "fundamental_cubic_matches_squared_coupling": sp.expand(
                2 * G * (-L * G * F) - dg.subs(g, G * G)
            ),
            "new_parameter_shift_changes_one_loop_only_at_order_three": sp.expand(
                h * (1 + h * h * F) - h
            ).coeff(h, 2),
            "no_new_first_order_shift_square": sp.expand(
                (G + h * h * (-L * G * F)) ** 2
            ).coeff(h, 2)
            - dg.subs(g, G * G),
        },
    }
