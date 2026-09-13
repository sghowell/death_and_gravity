"""Connected affine domain and full complex-neighborhood coefficient bounds."""

from functools import cache

import sympy as s

from . import family

OUTER_RADIUS = s.Rational(1, 16)
INNER_RADIUS = s.Rational(1, 32)
COMPLEX_EXPONENT = s.Rational(21, 32)
MAX_ORDER = 4
TUBE_BOUND = s.Rational(1, 10**2700)


@cache
def data():
    A = family.LOCALIZER
    exp_bound = s.factorial(8) * 2**8 / A**8
    delta_global = s.Rational(2**18, A)
    edge_gap = s.Rational(5, 48) ** family.N / 32
    delta_edge = 2**18 * s.factorial(4) / A**4
    f_pref = (
        family.LAMBDA * family.K0 * s.Rational(19, 16) ** 2
        + abs((family.GAMMA / 3 + family.CONTACT / 24) * family.K0)
        * s.Rational(17, 16) ** 4
    )
    j_pref = family.G * s.sqrt(family.K0) * s.Rational(17, 16) ** 2 / 2
    a, r, p, d, e, X = s.symbols("a R R_X delta_R delta_R_X X", real=True)
    qdiff = (p + e) ** 2 / (r + d) - p * p / r
    checks = {
        "outer_disk_real_square_minimum": s.Rational(13, 16) ** 2
        - s.Rational(1, 16) ** 2
        - COMPLEX_EXPONENT,
        "outer_switch_ratio": s.Rational(3, 16) / s.Rational(13, 16)
        - s.Rational(3, 13),
        "complex_base_R_deviation_product": s.Rational(17, 16)
        * s.Rational(3, 16)
        * s.Rational(4, 3)
        - s.Rational(17, 64),
        "full_inverse_ratio_difference": s.cancel(
            1 / (r + d) - 1 / r + d / (r * (r + d))
        ),
        "complete_Ia_square_difference": s.cancel(
            qdiff - ((2 * p * e + e * e) * r - p * p * d) / (r * (r + d))
        ),
        "second_order_full_product_norm": 2 * 97 * 65536 + 4 * 48**2 * 2048 - 31588352,
        "degree_four_Cauchy_factor": s.factorial(4) * 32**4 - 25165824,
        "uniform_lower_ratio_family": s.cancel(1 + (r - 1) / a - ((1 - 1 / a) + r / a)),
    }
    # Generic-R covariant metric map and affine mass factor.
    y = s.Symbol("positive_sqrt_R", positive=True)
    C = 1 / y
    D = (1 - C) / X
    checks["metric_norm_map"] = s.cancel(C + X * D - 1)
    checks["metric_inverse_map"] = s.cancel((1 - 1 / C) / X + D / C)
    checks["regular_null_metric_coefficient"] = s.cancel(
        D - (y * y - 1) / (X * y * (y + 1))
    )
    paff = s.Symbol("p", positive=True)
    tau = 3 * (2 * paff**3 - 1) / paff
    sigma = (8 * paff + 5) / (8 * paff**2)
    xi_t, xi_s = 1 - 1 / tau, -1 - 1 / sigma
    checks["full_affine_mass_determinant_factor"] = s.cancel(
        (1 + xi_t * tau) * (1 + xi_s * sigma) ** 3 + tau * sigma**3
    )
    schur = 1 - 3 * (r - 1) ** 2 / (2 * r)
    checks["temporal_schur_margin"] = s.cancel(
        schur - s.Rational(1, 4) - 3 * (r - s.Rational(1, 2)) * (2 - r) / (2 * r)
    )
    q, ru, H, Box, Z = s.symbols("q_regular R_u Href Box_u Z_u", real=True)
    alpha = 3 * (q + ru / 2) / 2 + 3 * (r - 1) * ru / (4 * r)
    cc = -(r - 1) / X
    dd = -(r - 1) / X**2 + 3 * (r - 1) * p / (2 * r * X)
    shift = 3 * H * (r - 1) + 3 * (q + ru / 2) / 2
    shifted = (r - 1) * (
        -3 * H + 3 * ru / (4 * r) + Box / X + (-1 / X**2 + 3 * p / (2 * r * X)) * Z
    )
    checks["new_q_recomputed_not_clock_reset"] = s.cancel(
        alpha - shift - cc * Box + dd * Z - shifted
    )
    checks["regular_q_cancels_from_retained_source"] = s.diff(shifted, q)
    aa, ap, bb, bp = s.symbols("a a_u b b_u", real=True)
    rg = 1 + aa * X**2 + bb * X**3
    rug = ap * X**2 + bp * X**3
    q4 = aa * ap / 3
    qjet = q4 * X**4
    ode = (
        s.diff(qjet, X)
        + (1 / (2 * X) - 3 * s.diff(rg, X) / (4 * rg)) * qjet
        - 3 * s.diff(rg, X) * rug / (4 * rg)
    )
    checks["new_regular_q_leading_integral_coefficient"] = (
        s.series(ode, X, 0, 4).removeO().expand().coeff(X, 3)
    )
    derivative_factors = {
        f"{i}_{j}": s.factorial(i) * s.factorial(j) * 32 ** (i + j)
        for i in range(5)
        for j in range(5 - i)
    }
    gates = {
        "unchanged_global_lower_via_nonnegative_delta": delta_global > 0,
        "global_low_X_branch_stays_below_six_fifths": s.Rational(9, 8) + delta_global
        < s.Rational(6, 5),
        "global_upper_edge_delta_below_half_literal_margin": delta_edge < edge_gap / 2,
        "all_added_R_values_tiny": delta_global < s.Rational(1, 10**414),
        "complex_bump_majorant_below_one": 2048
        * s.factorial(3)
        / (COMPLEX_EXPONENT * family.N) ** 3
        < 1,
        "complex_switch_error_below_one_sixteenth": 4 * s.Rational(3, 13) ** family.N
        < s.Rational(1, 16),
        "complex_scalar_denominator_bound": s.Rational(256, 255) ** 3
        < s.Rational(4, 3),
        "complex_old_R_modulus_between_half_and_three_halves": s.Rational(17, 64)
        < s.Rational(1, 2),
        "complex_new_R_stays_nonzero": 2048 * exp_bound < s.Rational(1, 16),
        "complex_old_R_X_Cauchy_bound": s.Rational(3, 2) * 32 == 48,
        "complex_delta_R_X_below_one": 65536 * exp_bound < 1,
        "complex_new_R_modulus_below_three_halves": s.Rational(17, 64)
        + s.Rational(1, 16)
        < s.Rational(1, 2),
        "inner_inverse_X_modulus_below_two": s.Rational(32, 27) < 2,
        "full_Ia_square_difference_below_two_power_25": 31588352 < 2**25,
        "full_Ia_A4_difference_below_two_power_26": 2**17 + s.Rational(7, 4) * 2**25
        < 2**26,
        "full_normalized_F_prefactor_below_one_e597": f_pref < 10**597,
        "full_normalized_J_prefactor_below_one_e397": j_pref < 10**397,
        "all_inner_coefficients_below_common_prefactor": 2**26 < 10**597,
        "every_mixed_Cauchy_factor_below_one_e8": all(
            v < 10**8 for v in derivative_factors.values()
        ),
        "full_every_mixed_four_jet_difference_below_one_e_minus_2700": 10**605
        * exp_bound
        < TUBE_BOUND,
        "affine_p_lower_from_original_R_lower": s.Rational(1, 2)
        > 4 * s.Rational(1, 3) ** 2,
        "affine_p_upper_from_original_R_upper": s.Rational(6, 5)
        < 4 * s.Rational(3, 5) ** 2,
        "affine_temporal_factor_negative": 2 * s.Rational(3, 5) ** 3 < 1,
    }
    return {
        "full_global_domain": "Every real u, -1/4096<X<6/5. The unchanged strict lower 1/2+1/(8*1024)<R_new and upper R_new<6/5 both hold. The sharper upper-edge estimate, not a tiny additive bound alone, retains the original affine domain.",
        "global_R_addition_bound": delta_global,
        "upper_edge_original_margin": edge_gap,
        "upper_edge_R_addition_bound": delta_edge,
        "tube": "Real |u|<=1 and 7/8<=X<=9/8. The estimate compares the full normalized coefficient functions F,R,A3,A4,A5 and J_H/sqrt(kappa0), not an arbitrary full solution.",
        "outer_joint_complex_radii": [OUTER_RADIUS, OUTER_RADIUS],
        "outer_real_X_squared_lower": COMPLEX_EXPONENT,
        "common_complete_coefficient_prefactor": s.Integer(10) ** 597,
        "inner_joint_complex_radii": [INNER_RADIUS, INNER_RADIUS],
        "all_mixed_derivative_factors": derivative_factors,
        "complete_actual_four_jet_majorant": s.Integer(10) ** 605 * exp_bound,
        "strict_four_jet_bound": TUBE_BOUND,
        "affine_construction": "Substitute the new full R into the generic S174 principal dictionary, the unique vacuum-regular q integral, the metric chart and the source-centered full rank60 mass update Xi=eta-D_old^-1. Keep the 56 complement equations and four projective gauges. Add the scalar/heavy action independently of the connection. Recompute B=3Href(R-1)+3(q+R_u/2)/2 and W=T-Bdu; the physical Proca source is independent of q after this exact centering. The affine background trace and coefficient maps need not equal their old values.",
        "regular_q_integral": "q=(3X/4)R(u,X)^(3/4) integral_0^1 sqrt(t) R_X(u,tX)R_u(u,tX)R(u,tX)^(-7/4)dt. Since R-1=X^2 a(u,X), this is regular at either sign of X and q=a*a_u*X^4/3+O(X^5). Its clock value is not reset.",
        "metric_map": "g_phys=R^-1/2 g_hat+(1-R^-1/2)du du/X with X_phys=X_hat, determinant ratio R^-3/2 and field-space Jacobian R^-9/2. All divided coefficients are regular at X=0, including nonzero null gradients.",
        "full_affine_mass_ratio": -tau * sigma**3,
        "reduced_source_after_new_centering": shifted,
        "scope": "This is a full connected algebraic parent and a quantitative coefficient tube. It is not a complete physical constraint count or stable cone on every background, an inverse estimate for a full quantum equation, or a nonlinear same-state bounce/UV/Regge theorem.",
        "checks": {key: s.cancel(value) for key, value in checks.items()},
        "gates": gates,
    }
