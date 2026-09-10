"""Finite-contact insertion as a full fixed-parameter one-loop variation."""

from functools import cache

import sympy as sp


@cache
def data():
    L, g = sp.symbols("quartic_L cubic_squared", real=True)
    h = sp.Symbol("external_heavy_inverse")
    t1, t2 = sp.symbols("first_internal_heavy_vertex second_internal_heavy_vertex")
    TA, TB, box, IR = sp.symbols(
        "first_integrated_heavy_triangle second_integrated_heavy_triangle integrated_heavy_box finite_light_bubble"
    )
    C = -L + g * h
    finite = C**2 * IR / 2 + C * (TA + TB) / 2 + box / 2
    derivative = sp.diff(finite, L)
    contact = sp.Symbol("once_fixed_finite_contact", real=True)
    x, z = sp.symbols("unit_Feynman_parameter diagnostic_channel")
    a = x * (1 - x)
    diagnostic = -L * sp.log(1 - a * z) / (16 * sp.pi**2)
    coeff = sp.integrate(sp.diff(diagnostic, z, 2).subs(z, 0) / 2, (x, 0, 1))
    ell, lp, cp = sp.symbols(
        "positive_logarithm_bound positive_quartic positive_cubic_over_mass",
        positive=True,
    )
    majorant = 3 * (4 * lp + 8 * cp * ell)
    return {
        "complete_integrated_one_loop_channel_without_constant_contact": finite,
        "fixed_g_and_M_quartic_variation": derivative,
        "once_fixed_contact_inserted_channel": contact * derivative,
        "complete_three_channel_derivative_majorant_without_loop_measure": (
            12 + 8 * ell
        )
        * lp,
        "diagnostic_nonzero_local_contact_loop_coefficient": contact * coeff,
        "scope": "Vary only L in the complete renormalized one-loop amplitude, with fixed g and M and all internal heavy propagators retained. The potential contact within that one-loop amplitude has zero second momentum derivative. The diagnostic pure-quartic bubble only demonstrates that a generic local contact need not remain momentum independent after insertion; it is not the actual model with g set to zero.",
        "checks": {
            "full_two_vertex_product_remainder_decomposition": sp.expand(
                (C + t1) * (C + t2) - C**2 - C * (t1 + t2) - t1 * t2
            ),
            "pointwise_quartic_variation_keeps_both_heavy_terms": sp.expand(
                sp.diff((C + t1) * (C + t2) - C**2, L) + t1 + t2
            ),
            "integrated_quartic_variation_keeps_both_triangles": sp.expand(
                derivative + C * IR + (TA + TB) / 2
            ),
            "heavy_box_has_no_quartic_variation_at_fixed_g_and_M": sp.diff(box, L),
            "bare_finite_contact_has_zero_second_momentum_derivative": sp.diff(
                -contact, z, 2
            ),
            "diagnostic_loop_contact_second_coefficient": sp.factor(
                coeff - L / (960 * sp.pi**2)
            ),
            "all_channel_derivative_majorant": sp.expand(
                majorant.subs(cp, lp / 3) - (12 + 8 * ell) * lp
            ),
            "fixed_logarithm_derivative_prefactor": 12 + 8 * 600 - 4812,
        },
    }
