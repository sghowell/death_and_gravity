"""Four-dimensional heat-kernel weight and full tensor normalization check."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes

from . import projectors as pro


@cache
def data():
    R2, Ric2, Riem2 = s.symbols("R2 Ric2 Riem2")
    scalar = R2 / 72 + (Riem2 - Ric2) / 180
    vector = 4 * scalar - R2 / 6 + Ric2 / 2 - Riem2 / 12
    proca = s.expand(vector - scalar)
    c = s.Rational(13, 120)
    e = -s.Rational(7, 40)
    r = s.Rational(1, 72)
    cov = c * (Riem2 - 2 * Ric2 + R2 / 3) + e * (Riem2 - 4 * Ric2 + R2) + r * R2
    x = s.symbols("h0:6", real=True)
    H = s.Matrix([[x[0], x[1], x[2]], [x[1], x[3], x[4]], [x[2], x[4], x[5]]])
    D = s.Matrix([[1, 2, -1], [2, -3, 1], [-1, 1, 4]])
    G = s.Matrix([[2, -1, 3], [-1, 1, 2], [3, 2, -2]]) / 7
    eps, eta = s.symbols("epsilon eta", real=True)
    sigma = s.Symbol("sigma", positive=True)
    Riem_quad = sigma**2 * s.trace(H * H)
    Ric_quad = sigma**2 * (s.trace(H) ** 2 + s.trace(H * H)) / 4
    R_quad = sigma**2 * s.trace(H) ** 2
    Weyl_quad = s.expand(Riem_quad - 2 * Ric_quad + R_quad / 3)
    source = eps * D + eta * G
    sub = {
        x[0]: source[0, 0],
        x[1]: source[0, 1],
        x[2]: source[0, 2],
        x[3]: source[1, 1],
        x[4]: source[1, 2],
        x[5]: source[2, 2],
    }
    Weyl_mixed = (
        s.expand(Weyl_quad.subs(sub, simultaneous=True)).coeff(eps, 1).coeff(eta, 1)
    )
    scalar_mixed = (
        s.expand(R_quad.subs(sub, simultaneous=True)).coeff(eps, 1).coeff(eta, 1)
    )
    checks = {
        "full_vector_minus_scalar_heat_weight": s.expand(proca - cov),
        "complete_Proca_Riemann_coefficient": proca.coeff(Riem2) + s.Rational(1, 15),
        "complete_Proca_Ricci_coefficient": proca.coeff(Ric2) - s.Rational(29, 60),
        "complete_Proca_scalar_coefficient": proca.coeff(R2) + s.Rational(1, 8),
        "full_Weyl_quadratic_form": s.expand(
            Weyl_quad - sigma**2 * (s.trace(H * H) - s.trace(H) ** 2 / 3) / 2
        ),
        "spin2_curvature_Hessian_factor": s.expand(
            Weyl_mixed - sigma**2 * (s.trace(D * G) - s.trace(D) * s.trace(G) / 3)
        ),
        "spin0_curvature_Hessian_factor": s.expand(
            scalar_mixed - 6 * sigma**2 * s.trace(D) * s.trace(G) / 3
        ),
        "UV_spin2_cut_weight": s.simplify(
            32
            * s.pi**2
            * s.limit(pro.above_threshold_density(2) / pro.S**2, pro.S, s.oo)
            - c
        ),
        "UV_spin0_cut_weight": s.simplify(
            32
            * s.pi**2
            * s.limit(pro.above_threshold_density(0) / pro.S**2, pro.S, s.oo)
            - 6 * r
        ),
        "actual_canonical_two_metric_factors": 4 / modes.KAPPA
        - 4 * s.Rational(1, 10) ** 800,
    }
    return {
        "heat_weight": "For the four-dimensional minimal Proca determinant, a4(vector)-a4(scalar)=-Riem^2/15+29Ric^2/60-R^2/8=13C^2/120-7Euler/40+R^2/72, modulo the separately tracked total-derivative coefficient. The Euler term is retained here.",
        "tensor_normalization": "For a transverse metric Fourier direction, the compact flat Weyl-squared Hessian is s^2 P2 and the R-squared Hessian is6s^2 P0. The full spatial stress cut accordingly has UV weights a4_C/(32pi^2) and6a4_R/(32pi^2).",
        "canonical": "For h=sqrt(kappa)gamma/2 the entire flat external-metric current dispersion density is multiplied by4/kappa. The spin0 projector is not asserted to be a canonically reduced propagating scalar or a full mixed constraint-reduced norm. This change of external-metric normalization does not make the response uniformly small at arbitrary spectral energy because it still grows as s^2.",
        "comparison_boundary": "This is a four-dimensional ultraviolet tensor-weight and normalization comparison, not a dimensionally continued finite local action. It does not redo S193, discard evanescent Euler contacts, change its finite coefficients or match S198's full curved contact/endpoint sector.",
        "physical_boundary": "The flat reference vacuum is a conditional external-metric Gaussian benchmark. No complete parent vacuum scattering amplitude, interacting loop remainder, gravity IR/Regge relation or P8 closure is asserted.",
        "checks": checks,
        "gates": {
            "nonzero_complete_spin2_weight": c > 0,
            "nonzero_complete_spin0_weight": r > 0,
            "Euler_weight_retained": e != 0,
            "same_fixed_canonical_parent": modes.MASS == 1000
            and modes.KAPPA == 10**800,
        },
    }
