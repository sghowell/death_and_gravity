"""Exact local F_alpha sum C^2 change and its fixed coordinate representation."""

from functools import cache

import sympy as s


@cache
def data():
    L, g, G, M, nu, F = s.symbols("L g G M nu F_alpha")
    hs = s.symbols("h_s h_t h_u")
    tree = -L + g * sum(hs)
    changes = {"L": -3 * L**2 * F, "g": -2 * L * g * F, "M": -g * F, "G": -L * G * F}
    amplitude = (
        -changes["L"]
        + changes["g"] * sum(hs)
        - g * changes["M"] * sum(H**2 for H in hs)
    )
    D = M - 2
    squared = (-L + g / (D - nu)) ** 2 + (-L + g / (D + nu)) ** 2 + (-L + g / M) ** 2
    b2 = s.diff(squared, nu, 2).subs(nu, 0) / 2
    btree = 2 * g / D**3
    relative = -2 * L + 3 * g / D
    return {
        "fixed_second_coordinate_reference": changes,
        "full_three_channel_C_squared": squared,
        "relative_b2_multiplier": relative,
        "MS_minus_old_b2": F * b2,
        "checks": {
            "literal_local_reference_amplitude": s.expand(
                amplitude - F * sum((-L + g * H) ** 2 for H in hs)
            ),
            "fundamental_G_reference_is_linear": 2 * G * changes["G"]
            - changes["g"].subs(g, G**2),
            "full_C_squared_forward_coefficient": s.factor(b2 / btree - relative),
            "tree_parameter_variation_forward_coefficient": s.factor(
                changes["g"] * s.diff(btree, g)
                + changes["M"] * s.diff(btree, M)
                - F * b2
            ),
            "potential_contact_cannot_tune_b2": s.diff(changes["L"], nu, 2),
            "t_channel_retained_but_constant": s.diff((-L + g / M) ** 2, nu, 2),
            "same_unshifted_tree_constant": s.diff(tree, L) + 1,
        },
        "scope": "The direct finite local amplitude term and this coordinate reorganization are the same contribution. No new free finite counterterm or extra external LSZ factor is used. For the actual g/(M-2)<L/2, the relative multiplier is negative and its absolute value is below 2L.",
    }
