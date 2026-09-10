"""One fixed local parent subtraction, not an adjustable b2 contact."""

from functools import cache

import sympy as sp


@cache
def data():
    L, g, M, F, G = sp.symbols(
        "quartic_L cubic_squared_g M paired_zero_bubble_F G", real=True
    )
    hs, ht, hu = sp.symbols("h_s h_t h_u", real=True)
    h = (hs, ht, hu)
    delta_L, delta_g, delta_M = 3 * L**2 * F, 2 * L * g * F, g * F
    ct = -delta_L + delta_g * sum(h) - g * delta_M * sum(v * v for v in h)
    subtract = F * sum((-L + g * v) ** 2 for v in h)
    s, t, u = sp.symbols("s t u", real=True)
    tree = -L + g * sum(1 / (M - z) for z in (s, t, u))
    actual_ct = (
        delta_L * sp.diff(tree, L)
        + delta_g * sp.diff(tree, g)
        + delta_M * sp.diff(tree, M)
    )
    q = sp.Symbol("finite_quartic_contact", real=True)
    finite = sp.Symbol("finite_change_in_paired_bubble_reference", real=True)
    forward_nonlocal_change = finite * sum(
        (-L + g / (M - z)) ** 2 for z in (s, 0, 4 - s)
    )
    nonlocal_b2 = sp.diff(forward_nonlocal_change, s, 2).subs(s, 2) / 2
    return {
        "paired_zero_bubble_reference": F,
        "delta_L_order_two_sector": delta_L,
        "delta_g_order_two_sector": delta_g,
        "delta_M_order_two_sector": delta_M,
        "delta_G_order_two_sector": L * G * F,
        "full_parent_counterterm_amplitude": ct,
        "subtraction": "F is the common regulated paired-insertion bubble at zero external momentum, including its 1/Q normalization. Subtract F sum C(z)^2 before regulator removal. It is not assigned a standalone finite integral.",
        "canonical_conversion": "The order-two sector is defined in this fixed local parent subtraction scheme. Its finite conversion to the S6.133 common MS interaction boundary and all other second-order shifts remain separate, unevaluated tasks.",
        "finite_nonlocal_reference_change_b2": nonlocal_b2,
        "checks": {
            "all_channels_implemented_by_local_parent_counterterms": sp.expand(
                ct + subtract
            ),
            "literal_full_tree_parameter_variation": sp.factor(
                actual_ct - ct.subs(dict(zip(h, [1 / (M - z) for z in (s, t, u)])))
            ),
            "fundamental_cubic_sector_increment": sp.factor(
                delta_g.subs(g, G**2) - 2 * G * (L * G * F)
            ),
            "same_constant_field_quartic_counterterm": sp.factor(
                delta_L
                - 3 * delta_g / M
                + 3 * g * delta_M / M**2
                - 3 * F * (L - g / M) ** 2
            ),
            "only_pure_contact_has_zero_b2": sp.diff(q, s, 2),
            "finite_nonlocal_reference_change_not_silently_zero": sp.factor(
                nonlocal_b2
                - finite * (-4 * L * g / (M - 2) ** 3 + 6 * g**2 / (M - 2) ** 4)
            ),
        },
    }
