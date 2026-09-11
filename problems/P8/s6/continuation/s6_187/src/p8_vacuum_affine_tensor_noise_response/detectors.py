"""Two exact boundary compatibilities and the leading noise detector pairing."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_metric_noise import noise
from p8_vacuum_affine_proca_gaussian.bridge import KAPPA

from . import energy


@cache
def data():
    t = s.Symbol("t", real=True)
    a = (1 + t * t) ** 2
    H = s.diff(a, t) / a
    k = s.Symbol("k", nonnegative=True)
    f = s.Function("homogeneous")(t)
    v = s.Function("advanced_test")(t)
    L = lambda x: s.diff(x, t, 2) + 3 * H * s.diff(x, t) + k * k * x / a**2
    boundary = a**3 * (f * s.diff(v, t) - s.diff(f, t) * v)
    derivative = s.diff(boundary, t)
    a0 = s.Rational(25, 16)
    V0, V1, I0, I1 = s.symbols("v_initial vprime_initial moment0 moment1", real=True)
    inverse = s.Matrix([I1 / a0**3, -I0 / a0**3])
    detector_norm = energy.data()["complete_adjoint_stress_smearing_norm_upper"]
    raw = noise.data()["complete_stress_variance_coefficient_upper"]
    h_bound = detector_norm**2 * raw / KAPPA
    gamma_bound = 4 * h_bound / KAPPA
    checks = {
        "actual_weighted_Green_identity": s.simplify(
            derivative - a**3 * (f * L(v) - v * L(f))
        ),
        "weighted_Wronskian_conserved": s.simplify(
            s.diff(a0**3 / a**3, t) + 3 * H * a0**3 / a**3
        ),
        "both_moment_to_initial_data_inversion": s.simplify(
            inverse.subs({I0: -(a0**3) * V1, I1: a0**3 * V0}) - s.Matrix([V0, V1])
        ),
        "canonical_tensor_noise_display": detector_norm * s.Rational(1, 10**375)
        - 5 * s.Rational(1, 10**372),
        "dimensionless_metric_noise_display": 2
        * 5
        * s.Rational(1, 10**372)
        / s.sqrt(KAPPA)
        - s.Rational(1, 10**771),
    }
    return {
        "fundamental_modes": "f0(t0)=1,f0'(t0)=0 and f1(t0)=0,f1'(t0)=1; W=f0 f1'-f0' f1=a0^3/a^3",
        "compatibility": "For each momentum and TT component, integral a^3 f0 q=0 and integral a^3 f1 q=0; equivalently the advanced v=Gadv q has v(t0)=v'(t0)=0",
        "advanced_initial_data": inverse,
        "exact_response_observable": "For q in this class, integral a^3 q h=T[v]/sqrt(kappa) for the source-only tree equation Lh=T_TT/sqrt(kappa). All homogeneous tensor contributions vanish in this pairing.",
        "leading_h_variance_upper_per_D_squared": h_bound,
        "leading_gamma_variance_upper_per_D_squared": gamma_bound,
        "leading_h_standard_deviation_display": "5e-372 D(q)",
        "leading_gamma_standard_deviation_display": "1e-771 D(q)",
        "initial_state_boundary": "No uncorrelated zero quantum metric initial state, sharp switch of the interaction or unrestricted finite-initial-time stress convolution is imposed. This is a restricted detector observable of the source-only leading component, not a complete metric state.",
        "finite_order_boundary": "Full tensor self-energy, fixed-profile response contacts, intrinsic/cross correlations at higher order and the finite-coupling remainder have not been bounded or discarded from a full-solution claim.",
        "checks": checks,
        "gates": {
            "actual_h_variance_below_display": h_bound
            < (5 * s.Rational(1, 10**372)) ** 2,
            "actual_gamma_variance_below_display": gamma_bound
            < s.Rational(1, 10**1542),
            "positive_initial_Wronskian": a0**3 > 0,
        },
    }
