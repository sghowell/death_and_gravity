"""Zero-soft momentum kernel, exact MS outer reference and finite local map."""

from functools import cache

import sympy as sp


@cache
def data():
    e = sp.Symbol("epsilon", real=True)
    m, v = sp.symbols("m v", positive=True)
    p = sp.Rational(3, 2) - e
    rho = v ** (1 - e) * (1 - 4 * m * m / v) ** p
    derivative = (
        8
        * p
        * v ** (-e)
        * (1 - 4 * m * m / v) ** (p - 2)
        * ((2 * p - 1) * 4 * m * m / v - 1)
    )
    B = (
        sp.sqrt(sp.pi)
        * sp.exp(2 * sp.EulerGamma * e)
        * 4 ** (-e)
        * (1 - 2 * e / 3)
        * (1 - 4 * e)
        * sp.gamma(1 + e)
        * sp.gamma(1 + 2 * e)
        / ((1 - e) * (1 + 2 * e) * sp.gamma(sp.Rational(1, 2) + e))
    )
    L = sp.log(B)
    first = sp.simplify(sp.diff(L, e).subs(e, 0))
    second = sp.simplify(sp.diff(L, e, 2).subs(e, 0))
    hardfinite = sp.simplify(-3 * (first**2 + second))
    ell = sp.Symbol("log_m_squared", real=True)
    A = sp.exp(sp.EulerGamma * e) * sp.gamma(1 + e)
    soft_normal = sp.exp(ell * e) * ((12 - 32 * e + 16 * e * e) * A * A - 12 * A)
    softfinite = sp.simplify(sp.diff(soft_normal, e, 2).subs(e, 0) / 2)
    k = sp.Symbol("k", positive=True)
    triangle = (v * sp.log(v) - v + 1) / (v - 1) ** 2
    delta = (
        1 / v * ((1 - 1 / v) ** (-2) - 1) * (sp.log(v) - 1)
        + 1 / v**2 / (1 - 1 / v) ** 2
    )
    regulated = sp.gamma(e) * (
        1 / (v - 1) - (v ** (1 - e) - 1) / ((1 - e) * (v - 1) ** 2)
    )
    # Evaluate its removable singularity via numerator derivative, retaining Gamma(e).
    tfinite = sp.diff(regulated / sp.gamma(e), e).subs(e, 0)
    Lc, g, M, F, s, t, u = sp.symbols("L g M F s t u")
    C = lambda z: -Lc + g / (M - z)
    tree = -Lc + g * sum(1 / (M - z) for z in (s, t, u))
    variation = sp.diff(tree, Lc) * (-3 * Lc * F) + sp.diff(tree, g) * (-g * F)
    w = sp.Symbol("w")
    b2 = (
        sp.diff(
            (-F * sum(C(z) for z in (s, t, u))).subs({s: 2 + w, t: 0, u: 2 - w}), w, 2
        ).subs(w, 0)
        / 2
    )
    return {
        "zero_soft_kernel": "Gamma_0,MS(q^2)=Y d^2 f_MS(-q^2)/dm^2 at fixed mu, then mu=m. Gamma_0(0)=-64 N Y^2/Q.",
        "spectral_density": "nu(v)=12 C Y (2T/v-1)/sqrt(1-T/v), C=2NY/Q,T=4m^2; signed, not a positive spectral measure.",
        "kernel_majorant": "|Gamma_0(q^2)| <= (NY^2/Q)[136+48 log(1+q^2)].",
        "dimensionless_hard_reference": -6 * B / e**2,
        "dimensionless_soft_forest": soft_normal / e**2,
        "leading_finite_reference": sp.simplify(softfinite + hardfinite),
        "leading_double_pole": -6,
        "leading_simple_pole": 2,
        "exact_4D_triangle": triangle,
        "finite_correction": "delta F/(CY/Q)=-12 integral_0^1 (2z-1)/(z sqrt(1-z))*[h(rz)(-log(rz)-1)+rz A(rz)] dz; h=A-1,A=(1-k)^(-2),r=1/T.",
        "correction_upper": "|delta F|/(CY/Q) <=r[300+100 log T], 0<r<=1/16.",
        "common_MS_minus_subtracted_amplitude": -F * sum(C(z) for z in (s, t, u)),
        "matched_star_shifts": {
            "L": -3 * Lc * F,
            "g": -g * F,
            "M": 0,
            "G_relative": -F / 2,
        },
        "b2_conversion": b2,
        "checks": {
            "fixed_mu_mass_derivative_spectral_density": sp.simplify(
                sp.diff(rho, m, 2) - derivative
            ),
            "hard_normalization": sp.simplify(B.subs(e, 0) - 1),
            "hard_first_log_jet": first + sp.Rational(17, 3),
            "hard_second_log_jet": second + sp.Rational(103, 9) - sp.pi**2 / 3,
            "hard_finite": hardfinite + 62 + sp.pi**2,
            "soft_double_pole_cancellation": soft_normal.subs(e, 0),
            "soft_simple_pole": sp.simplify(sp.diff(soft_normal, e).subs(e, 0) + 32),
            "soft_finite": softfinite - 16 - sp.pi**2 + 32 * ell,
            "full_finite": sp.simplify(softfinite + hardfinite + 46 + 32 * ell),
            "full_simple_pole": 34 - 32 - 2,
            "triangle_removable_regulator_limit": sp.simplify(tfinite - triangle),
            "exact_finite_mass_ratio_difference": sp.simplify(
                triangle - (sp.log(v) - 1) / v - delta
            ),
            "h_below_three_k_polynomial": sp.factor(
                3 * k
                - ((1 - k) ** (-2) - 1)
                - k * (1 - 5 * k + 3 * k * k) / (1 - k) ** 2
            ),
            "A_below_two_on_domain": 2
            - sp.Rational(16, 15) ** 2
            - sp.Rational(194, 225),
            "finite_remainder_conservative_constant": 12 * (3 * 4 + 5 * 2) - 264,
            "finite_remainder_log_coefficient": 12 * 6 - 72,
            "conversion_literal_tree_derivative": sp.factor(
                variation + F * sum(C(z) for z in (s, t, u))
            ),
            "b2_relative_conversion": sp.factor(b2 + 2 * g * F / (M - 2) ** 3),
            "zero_soft_MS_reference": -32 * 2 + 64,
        },
    }
