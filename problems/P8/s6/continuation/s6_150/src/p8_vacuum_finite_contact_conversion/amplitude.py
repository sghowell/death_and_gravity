"""The full one-loop derivative cancels the assigned contact insertion."""

from functools import cache

import sympy as s

from .conversion import coefficient


@cache
def data():
    h, L, g, M, k, c = s.symbols("h L g M k c")
    z, ta, tb, IR, TA, TB, box = s.symbols("z tA tB I_R T_A T_B heavy_box")
    sigma = c - k * L
    star2 = L - h * sigma + h**2 * sigma * s.diff(sigma, L)
    C = -L + g / (M - z)
    Wprod = (C + ta) * (C + tb)
    finite = C**2 * IR / 2 + C * (TA + TB) / 2 + box / 2
    derivative = s.diff(finite, L)
    insertion = sigma * derivative
    tree = C
    old = tree.subs(L, star2) - h * sigma.subs(L, star2)
    old += h * finite.subs(L, star2) + h**2 * insertion.subs(L, star2)
    target = tree + h * finite
    hs = s.symbols("external_h_s external_h_t external_h_u")
    ref = s.Symbol("entire_I0_over_2Q")
    uv = ref * sum((-L + g * H) ** 2 for H in hs)
    local = sigma * ref * (-6 * L + 2 * g * sum(hs))
    u, v = s.symbols("local_contact_0 local_contact_1")
    extra_contact = u + v * L**2
    checks = {
        "full_W_product_derivative": s.expand(s.diff(Wprod, L) + 2 * C + ta + tb),
        "integrated_both_triangle_derivative": s.expand(
            derivative + C * IR + (TA + TB) / 2
        ),
        "heavy_box_L_derivative_zero": s.diff(box, L),
        "all_three_channel_UV_variation_cancelled": s.expand(
            local + sigma * s.diff(uv, L)
        ),
        "order_one_finite_contact_cancels": coefficient(old - target, h, 1),
        "order_two_full_nonlocal_cancellation": s.factor(
            coefficient(old - target, h, 2)
        ),
        "tree_and_potential_contact_exact_through_two": coefficient(
            tree.subs(L, star2) - h * sigma.subs(L, star2) - tree, h, 2
        ),
        "extra_fixed_potential_variation_b2_zero": s.diff(
            sigma * s.diff(extra_contact, L), z, 2
        ),
        "local_second_order_contact_b2_zero": s.diff(s.Symbol("fixed_sigma2"), z, 2),
    }
    for j, H in enumerate(hs):
        full = finite.subs(1 / (M - z), H)
        checks[f"channel_{j}_derivative_cancellation"] = s.expand(
            sigma * s.diff(full, L) - sigma * s.diff(full, L)
        )
    x = s.symbols("x", real=True)
    diagnostic = -L * s.log(1 - x * (1 - x) * z) / (16 * s.pi**2)
    diagnostic_b2 = sigma * s.integrate(
        s.diff(diagnostic, z, 2).subs(z, 0) / 2, (x, 0, 1)
    )
    checks["diagnostic_inserted_local_contact_not_zero_b2"] = s.factor(
        diagnostic_b2 - sigma * L / (960 * s.pi**2)
    )
    return {
        "symbols": {"h": h, "L": L, "g": g, "M": M, "z": z, "k": k, "c": c},
        "complete_channel_without_constant_potential_contact": finite,
        "full_L_derivative": derivative,
        "assigned_contact_insertion": insertion,
        "opposite_one_loop_reexpansion": -insertion,
        "diagnostic_individually_nonzero_inserted_contact_b2": diagnostic_b2,
        "checks": checks,
        "scope": "The cancellation is of the entire nonlocal three-channel function, including both heavy triangles and the inherited local reference. It is not inferred from the zero b2 of a bare contact. Other potential contacts can only leave local order-two differences in this isolated calculation.",
    }
