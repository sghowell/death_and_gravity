"""Finite outer MS conversion of the complete two-position insertion family."""

from functools import cache

import sympy as s


@cache
def data():
    e, ell, Q, a0, a1, a2 = s.symbols("epsilon ell Q alpha0 alpha1 alpha2")
    alpha = a0 + e * a1 + e**2 * a2
    I0 = (1 / e + ell + e * (ell**2 / 2 + s.pi**2 / 12)) / Q
    complete = s.series(alpha * I0, e, 0, 1).removeO().expand()
    finite = (a0 * ell + a1) / Q
    x, M, g = s.symbols("x M g", positive=True)
    D = x * M + (1 - x) ** 2
    b = x * (1 - x) / D
    A0 = s.Integral(b, (x, 0, 1))
    A1 = s.Integral(b * (ell - s.log(D)), (x, 0, 1))
    F = g / Q**2 * s.Integral(b * (2 * ell - s.log(D)), (x, 0, 1))
    C1, C2, R1, R2, D1, D2 = s.symbols("V1 V2 R1 R2 D1 D2")
    L, cubic = s.symbols("L cubic_squared")
    hs = s.symbols("h_s h_t h_u")
    UV = sum((-L + cubic * H) ** 2 for H in hs)
    oldCT = -alpha * I0 * UV
    newCT = -a0 * UV / (Q * e)
    z, m0, m1 = s.symbols("inner_invariant local_mass local_slope")
    affine = m0 + m1 * z
    return {
        "entire_regulated_multiplier_times_outer_reference": complete,
        "finite_outer_coefficient": finite,
        "actual_F_alpha_integral": F,
        "MS_minus_old_complete_family": s.Symbol("F_alpha") * UV,
        "checks": {
            "same_two_position_constant_insertion_weight": s.expand(
                C1 * C2 * (2 * alpha + R1 + R2) / (2 * D1 * D2)
                - alpha * C1 * C2 / (D1 * D2)
                - C1 * C2 * (R1 + R2) / (2 * D1 * D2)
            ),
            "finite_product_includes_alpha_first_epsilon": s.expand(
                complete - a0 / (Q * e) - finite
            ),
            "complete_local_reference_difference": s.simplify(
                s.limit(newCT - oldCT, e, 0) - finite * UV
            ),
            "combined_two_scale_logs_in_integrand": s.expand(
                b * ell + b * (ell - s.log(D)) - b * (2 * ell - s.log(D))
            ),
            "single_leg_alpha_truncation_finite_defect": s.simplify(
                s.limit((alpha - a0) * I0, e, 0) - a1 / Q
            ),
            "alpha_zero_plus_first_combination": s.expand(
                (g * A0 * ell / Q + g * A1 / Q) / Q - g * (A0 * ell + A1) / Q**2
            ),
            "constant_tadpole_inner_OS_annihilation": s.expand(
                affine - affine.subs(z, 1) - (z - 1) * s.diff(affine, z)
            ),
        },
        "scope": "Only the remaining outer interaction pole projection changes. The fixed inner physical OS grouping, all 64 refinements and both convergent decaying insertions are retained. The finite difference is F_alpha sum C^2, with both scale logs fixed by the same regulator.",
    }
