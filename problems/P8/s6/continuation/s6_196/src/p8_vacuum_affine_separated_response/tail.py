"""Quantitative all-momentum tail of the actual centered stress vector."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_metric_noise import oscillatory, reference
from p8_vacuum_affine_spatial_current.bounds import require_band
from p8_vector_state.comparison import AMAX

DISPLAY = s.Integer(10) ** 52


def variance_tail_display(K):
    return DISPLAY / require_band(K)


@cache
def constants():
    A, m = AMAX, reference.MASS
    ref = oscillatory.REF_BOUND**2 * A**4
    radial10 = 64 * A**10 / 63
    low = 36 * reference.PAIR_ERROR**2 * 8 * radial10
    high = 36 * reference.PAIR_ERROR**2 * 4 / (100 * m**7)
    return {
        "reference_tail_variance_numerator": ref,
        "tenth_radial_half_band_numerator_pi_lower_bound": radial10,
        "low_external_momentum_error_tail_numerator": low,
        "high_external_momentum_error_tail_numerator": high,
        "high_external_weight_bound_numerator": s.Rational(4, 1000) + 2 / (A * m),
        "uniform_total_tail_numerator_before_rounding": ref
        + low / s.Integer(1000) ** 6
        + high,
    }


@cache
def data():
    K, x, p = s.symbols("K x p", positive=True)
    A, m = AMAX, reference.MASS
    nu, mu = s.symbols("nu mu", positive=True)
    c = constants()
    L = 1 + p / (A * m)
    checks = {
        "fourth_radial_infinite_tail": s.integrate(x**-2, (x, K, s.oo)) - 1 / K,
        "tenth_radial_infinite_half_band_tail": s.integrate(x**-8, (x, K / 2, s.oo))
        - 128 / (7 * K**7),
        "reference_AM_GM": s.expand((nu + mu) ** 2 - 4 * nu * mu - (nu - mu) ** 2),
        "complete_remainder_square_split": s.factor(
            2 * (nu**-12 + mu**-12) - (nu**-6 + mu**-6) ** 2 - (nu**-6 - mu**-6) ** 2
        ),
        "external_linear_to_H1_weight": s.expand(
            2 * (1 + p * p)
            - L
            - 2 * (p - 1 / (4 * A * m)) ** 2
            - (1 - 1 / (8 * A * A * m * m))
        ),
        "actual_fixed_mass": m - 1000,
        "canonical_partition_value": variance_tail_display(10**16) - 10**36,
    }
    return {
        "actual_two_particle_tail": "tau_f(K)^2=||(1-P_K^(2)) T_centered[f]Omega||^2. The same actual pure all-order state and all nine polarization pairs are used; P_K^(2) retains BOTH created momenta |k|,|l|<=K.",
        "reference": "The tail union max(|k|,|l|)>K is bounded by J4_tail(K)/2 uniformly in P=k+l. With J4_tail<=Amax^4/(2pi^2 K), the complete reference variance numerator is below1e51.",
        "low_external": "For |P|<=K/2 the tail union forces both legs>K/2. The complete remainder integral is <=4L(P)J10_tail(K/2), with J10_tail<100/K^7 and L<2(1+|P|^2). Its variance numerator is below1e35/K^7.",
        "high_external": "For |P|>K/2 keep the full internal remainder integral4L(P)J10 and use L(P)/(1+|P|^2)<=4/K^2+2/(Amax m K)<1/(100K). Its variance coefficient is below1e10/K. No external momentum cutoff is imposed.",
        "result": "For K>=1000, tau_f(K)^2<1e52 N[f]^2/K for nonzero f, and it is zero when f=0. N is exactly the S186 three-time/one-spatial-derivative L2 tensor norm.",
        "constants": c,
        "checks": checks,
        "gates": {
            "reference_complete_pair_bound": c["reference_tail_variance_numerator"]
            < 10**51,
            "exact_half_band_tenth_radial_display": c[
                "tenth_radial_half_band_numerator_pi_lower_bound"
            ]
            < 100,
            "low_external_remainder_display": c[
                "low_external_momentum_error_tail_numerator"
            ]
            < 10**35,
            "high_external_H1_weight": c["high_external_weight_bound_numerator"]
            < s.Rational(1, 100),
            "high_external_remainder_display": c[
                "high_external_momentum_error_tail_numerator"
            ]
            < 10**10,
            "complete_tail_after_K_at_least_mass": c[
                "uniform_total_tail_numerator_before_rounding"
            ]
            < DISPLAY,
            "positive_external_weight_floor": 1 - 1 / (8 * A * A * m * m) > 0,
        },
    }
