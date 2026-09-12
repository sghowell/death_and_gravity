"""All-spatial-momentum coefficient and spectral-tail bounds."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_flat_tensor_cut import dispersion, projectors

MASS = projectors.MASS


@cache
def data():
    m = MASS
    sigma, q = s.symbols("sigma q", positive=True)
    total = dispersion.moment(2, 3) + dispersion.moment(0, 3)
    spectral_tail = (s.Rational(1, 128) + s.Rational(1, 384)) / s.pi**2
    checks = {
        "complete_sine_tensor_density_weight": 2 * s.Rational(1, 128)
        + s.Rational(1, 384)
        - s.Rational(7, 384),
        "absolute_sine_tail_antiderivative": s.diff(-2 / s.sqrt(sigma + q), sigma)
        - (sigma + q) ** (-s.Rational(3, 2)),
        "complete_absolute_sine_radial_constant": 2 * s.Rational(7, 384)
        - s.Rational(7, 192),
        "canonical_full_flat_nonlocal_display": 4 * (2 * 10**10) / modes.KAPPA
        - 8 * s.Rational(1, 10) ** 790,
        "canonical_full_flat_nonlocal_tail": 4 * 10**14 / modes.KAPPA
        - 4 * s.Rational(1, 10) ** 786,
        "complete_first_moment_sum": s.simplify(total - 3 / (2560 * s.pi**2 * m * m)),
        "coefficient_Cauchy_bound": s.simplify(
            s.Rational(32, 3) * total - 1 / (80 * s.pi**2 * m * m)
        ),
        "complete_spectral_tail_constant": s.simplify(
            spectral_tail - 1 / (96 * s.pi**2)
        ),
        "coefficient_Cauchy_tail": s.simplify(
            s.Rational(32, 3) * spectral_tail - 1 / (9 * s.pi**2)
        ),
        "canonical_standard_norm_display": 4 * 10**10 / modes.KAPPA
        - 4 * s.Rational(1, 10) ** 790,
        "canonical_standard_norm_tail_display": 4 * 10**17 / modes.KAPPA
        - 4 * s.Rational(1, 10) ** 783,
    }
    return {
        "analytic_disc": "For fixed real p, q=|p|^2, use |z+q|<=m^2. Then |s-z|>=3s/4 for s>=4m^2, |z|<=q+m^2 and ||z I+p p^T||op<=2(q+m^2).",
        "full_tensor_norm": "On complex symmetric spatial matrices with Frobenius norm, ||Q0||<=||A||op^2 and ||Q2||<=2||A||op^2, using the rank-one trace term as well as A D A^T. Thus ||F||<=(32/3)(I2+I0)(q+m^2)^3 on the disc.",
        "coefficient_bound": "Cauchy gives ||A_r(p)||<=(q+m^2)^3/[80pi^2 m^(2r+2)] for r=0,1,2, at every real spatial momentum.",
        "tail_bound": "With a computational squared-spectral-mass limit Lambda>=4m^2, ||A_r-A_r,Lambda||<=(q+m^2)^3/[9pi^2 m^(2r) Lambda]. This bounds actual improper tails, not an assumed physical cutoff.",
        "weighted_norm": "Let Y(Gamma)^2=sum_(r=0)^2 m^(-4r)||(m^2-Delta)^3 partial_t^(2r) Gamma||L2^2. The conversion bilinear is bounded by sqrt(3)/(80pi^2 m^2)||D||L2 Y, with tail sqrt(3)/(9pi^2 Lambda)||D||L2 Y.",
        "standard_norm": "For m=1000, let Z(Gamma)^2=sum_(r=0)^2||(1-Delta)^3 partial_t^(2r) Gamma||L2^2. Since m>=1, the full conversion is below1e10||D||L2 Z and its spectral tail below1e17||D||L2 Z/Lambda. The source norm includes six spatial and four time derivatives.",
        "canonical": "Both external-metric chain factors give4e-790||D||L2 Z and4e-783||D||L2 Z/Lambda. These derivative-losing flat-benchmark bounds are not a fully reduced mixed norm, same-space inverse or contraction.",
        "complete_flat_bulk": "The actual sixth-derivative spectral sine bulk has matrix majorant integral (2rho2+rho0)(1+q/s)^2/(s+q)^(7/2) ds <=7/[192pi^2 sqrt(4m^2+q)]. On the unit slab its bilinear is below1e-5||D||L2||partial_t^6 Gamma||L2 for m=1000, with spectral tail below||D||L2||partial_t^6 Gamma||L2/(100sqrt(Lambda)).",
        "full_nonlocal_response": "Let U(Gamma)^2=Z(Gamma)^2+||partial_t^6 Gamma||L2^2. The complete specified flat covariant nonlocal representative has weak bound2e10||D||L2 U, and full spectral tail1e14||D||L2 U/sqrt(Lambda) for Lambda>=4m^2. Canonical displays are8e-790 and4e-786/sqrt(Lambda). The physical local polynomial and full curved response are not included.",
        "checks": checks,
        "gates": {
            "actual_mass_and_normalization": modes.MASS == 1000
            and modes.KAPPA == 10**800,
            "ordinary_conversion_display": s.sqrt(3) * 1000**4 / (80 * s.pi**2)
            < 10**10,
            "ordinary_tail_display": s.sqrt(3) * 1000**6 / (9 * s.pi**2) < 10**17,
            "positive_first_moments": total > 0,
            "flat_bulk_unit_slab_display": 7 / (384 * s.pi**2 * 1000)
            < s.Rational(1, 10**5),
            "flat_bulk_tail_display": 7 / (192 * s.pi**2) < s.Rational(1, 100),
            "full_flat_nonlocal_display": 10**10 + s.Rational(1, 10**5) < 2 * 10**10,
            "full_flat_nonlocal_tail_display": s.Rational(10**17, 2000)
            + s.Rational(1, 100)
            < 10**14,
        },
    }
