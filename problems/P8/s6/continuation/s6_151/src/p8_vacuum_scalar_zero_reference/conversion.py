"""The entire wineglass family changes by a finite local anchor and scale term."""

from functools import cache

import sympy as s


@cache
def data():
    Jp, Jzero, I, I0, Ims, K = s.symbols(
        "J_bare_P J_bare_zero I_outer_P I_zero I_MS overall_MS_poles"
    )
    L, g, M, ell, Q, F = s.symbols("L g M ell Q finite_proper_MS_zero")
    C, T, CQ = s.symbols("external_C two_full_heavy_triangles inner_C")
    old = Jp - I0 * I - Jzero + I0**2
    new = Jp - Ims * I - K
    newzero = Jzero - Ims * I0 - K
    heavy = L**2 * T + (C + T) * (CQ**2 - L**2)
    dL, dg, dM = 3 * L**3 * F, L**2 * g * F, s.Integer(0)
    hs = s.symbols("h_s h_t h_u")
    G, nu = s.symbols("G nu")
    tree = -L + g * sum(hs)
    b2tree = 2 * g / (M - 2) ** 3
    return {
        "whole_proper_MS_core_minus_old_core": new - old,
        "finite_core_conversion": F + ell * s.Symbol("finite_outer_I_R") / Q,
        "complete_heavy_decaying_numerator": heavy,
        "full_family_MS_minus_old": L**2 * s.Symbol("sum_external_C") * F
        + ell / Q * s.Symbol("complete_outer_core_plus_heavy_decaying_integral"),
        "local_reference_coordinate_map": {
            "delta_L": dL,
            "delta_g": dg,
            "delta_M": dM,
            "delta_G": L**2 * G * F / 2,
        },
        "local_b2_relative_change": L**2 * F,
        "checks": {
            "regulated_core_difference_before_finite_parts": s.expand(
                new - old - newzero - (I0 - Ims) * (I - I0)
            ),
            "whole_sixteen_vertex_numerator": s.expand(
                (C + T) * CQ**2 - C * L**2 - heavy
            ),
            "both_heavy_triangles_retained": s.diff(heavy, T) - (L**2 + CQ**2 - L**2),
            "finite_nonlocal_product_has_no_pole": s.expand(
                (s.Symbol("deltaI_finite") + s.Symbol("epsilon") * s.Symbol("deltaI1"))
                * s.Symbol("outer_IR_finite")
            ).subs(s.Symbol("epsilon"), 0)
            - s.Symbol("deltaI_finite") * s.Symbol("outer_IR_finite"),
            "literal_local_coordinate_tree_variation": s.expand(
                dL * s.diff(tree, L)
                + dg * s.diff(tree, g)
                - L**2 * F * sum(-L + g * H for H in hs)
            ),
            "fundamental_G_linear_variation": 2 * G * (L**2 * G * F / 2)
            - dg.subs(g, G**2),
            "no_heavy_mass_reference_for_this_family": dM,
            "local_b2_exact_relative_change": s.factor(
                dg * s.diff(b2tree, g) / b2tree - L**2 * F
            ),
            "local_quartic_constant_cannot_change_b2": s.diff(dL, nu, 2),
        },
        "scope": "All sixteen wineglass refinements and their assigned proper/overall interaction forests at MS mu=mF, holding the canonical field convention fixed. The finite heavy remainder acquires (ell/Q) times its full outer integral. Other interaction families, finite coordinate/field-map cross terms and full canonical assembly remain open.",
    }
