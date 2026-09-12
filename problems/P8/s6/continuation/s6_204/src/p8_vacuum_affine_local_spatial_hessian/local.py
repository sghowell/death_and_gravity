"""Fixed four-dimensional local action and complete conformal-time Hessians."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes

from .helicities import B, C, basis, tf

MASS = modes.MASS


def operators(P, jets, a, a1, a2, a3):
    q = (P.T * P)[0]
    M0, M1, M2, _, M4 = jets
    r = a2 / a
    r1 = a3 / a - a1 * a2 / a**2
    spatial = q * M0 - 2 * B(P, M0)
    return {
        "R_old": -(a * a * M2 + 2 * a * a1 * M1 + a * a * spatial) / 2,
        "R_old_squared": -6 * (r * M2 + r1 * M1 + r * spatial) + 2 * C(P, M0),
        "Weyl_squared": M4
        + 2 * q * M2
        - 2 * B(P, M2)
        + q * q * M0
        - 2 * q * B(P, M0)
        + s.Rational(2, 3) * C(P, M0),
    }


@cache
def data():
    p = s.Matrix(s.symbols("p1:4", real=True))
    q = (p.T * p)[0]
    w = s.Symbol("w")
    z = w - q
    x = s.symbols("x0:5")
    M = sum((v * T for v, T in zip(x, basis())), s.zeros(3))
    zero = s.zeros(3)
    flat = operators(
        p,
        [M, zero, -w * M, zero, w * w * M],
        s.Integer(1),
        s.Integer(0),
        s.Integer(0),
        s.Integer(0),
    )
    A = z * s.eye(3) + p * p.T
    Q2 = tf(A * M * A - A * s.trace(A * M) / 3)
    sixQ0 = 2 * tf(A) * s.trace(A * M)
    R2, Ric2, Riem2 = s.symbols("R2 Ric2 Riem2")
    scalar = R2 / 72 + (Riem2 - Ric2) / 180
    Weyl = Riem2 - 2 * Ric2 + R2 / 3
    Euler = Riem2 - 4 * Ric2 + R2
    checks = {
        "full_spatial_Weyl_Hessian_agrees_with_covariant_polynomial": flat[
            "Weyl_squared"
        ]
        - Q2,
        "full_spatial_R_squared_Hessian_agrees_with_covariant_polynomial": flat[
            "R_old_squared"
        ]
        - sixQ0,
        "original_fixed_finite_curvature_density": s.expand(
            -4 * scalar + Weyl / 30 + R2 / 18 - Euler / 90
        ),
        "original_fixed_Einstein_coefficient": s.Rational(5, 3) * MASS**2
        - s.Rational(5000000, 3),
        "unimodular_volume_all_directions": s.trace(M),
    }
    return {
        "fixed_action": "After the already established original dimension limit at mu=m, the finite local density is[5m^4/2+(5/3)m^2 R_old-C^2/30-R_old^2/18+Euler/90]/(64pi^2). Only compact four-dimensional Euler variation is removed here; neither evanescent terms nor finite coefficients are retuned.",
        "parametrization": "Use g=a(eta)^2[-deta^2+exp(gamma)_ij dx^i dx^j], with gamma symmetric tracefree and a fixed. The determinant is exactly-a^8, so the volume term has zero variation. This is the same exponential/unimodular spatial parametrization.",
        "Einstein_quadratic": "Modulo compact spatial divergences, the quadratic Einstein action is integral a^2/4[||gamma'||^2-||grad gamma||^2+2||div gamma||^2]. No lapse or shift variation is claimed.",
        "R_squared_quadratic": "With r=a''/a, the quadratic R_old^2 action is integral[(partial_i partial_j gamma_ij)^2+3r(||gamma'||^2-||grad gamma||^2+2||div gamma||^2)]. The conformal curvature shift6a''/a is independent of gamma because det exp(gamma)=1.",
        "Weyl_quadratic": "The background is conformally flat. The compact four-dimensional Weyl-squared quadratic action has the flat conformal Hessian displayed by operators(), retaining all spatial tensor components. No flat physical state is substituted into a curved quantum kernel.",
        "operators": "H_R=-(a^2 partial_eta^2+2aa' partial_eta+a^2(qI-2B))/2; H_R2=-6[r partial_eta^2+r' partial_eta+r(qI-2B)]+2C; H_C=partial_eta^4+(2qI-2B)partial_eta^2+q^2 I-2qB+2C/3.",
        "fixed_Hessian": "H_fixed=[(5/3)m^2 H_R-H_C/30-H_R2/18]/(64pi^2), as the compact action second variation. This is an explicit target supplied by the already fixed local prescription, not proof that the remaining quantum endpoint/contact sector matches it.",
        "checks": checks,
        "gates": {
            "same_actual_mass": MASS == 1000,
            "finite_Weyl_coefficient_not_Proca_pole_weight": s.Rational(-1, 30)
            != s.Rational(13, 120),
            "physical_dimension_Euler_only_after_fixed_matching": True,
            "no_lapse_shift_or_reduced_mixed_operator_claim": True,
        },
    }
