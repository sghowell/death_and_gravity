"""Both clock-normalized and canonical source-force covariance conventions."""

from functools import cache

import sympy as s
from p8_vacuum_affine_retarded_energy import estimates as previous

from . import derivatives, modes

SMALL_DELTA = s.Rational(1, 10**14)


@cache
def data():
    m, K = modes.MASS, modes.KAPPA
    d = s.Symbol("delta", positive=True)
    coeff = 108 * m * 2048**2 + 144 * derivatives.GRAD**2 / m
    normalized = coeff * d * d / K
    canonical = coeff * d * d
    mean = s.Integer(4) * 10**10 * d * d
    field, test = s.symbols(
        "canonical_perturbation_norm canonical_test_norm", nonnegative=True
    )
    raw_psi, raw_eta = field / s.sqrt(K), test / s.sqrt(K)
    checks = {
        "complete_covariance_normalization": normalized
        - (m * m / K) * (108 * 2048**2 / m + 144 * derivatives.GRAD**2 / m**3) * d * d,
        "canonical_test_map_cancels_normalized_gravity_smallness": K
        * K
        * (coeff * d * d / K)
        * (1 / K)
        - canonical,
        "canonical_mean_force_pairing_dictionary": K
        * (4 * 10**10 * d * d)
        * raw_psi
        * raw_eta
        - mean * field * test,
        "source_sector_same_previous_mean_bound": previous.data()["kappa"] - K,
        "reference_noise_zero": normalized.subs(d, 0),
        "reference_canonical_noise_zero": canonical.subs(d, 0),
        "canonical_small_ball_standard_deviation_display": 2 * 10**8 * SMALL_DELTA
        - 2 * s.Rational(1, 10**6),
        "canonical_small_ball_mean_relative_display": 4 * 10**10 * SMALL_DELTA**2
        - 4 * s.Rational(1, 10**18),
    }
    return {
        "noise_identity": "(F_eta-mean)/kappa=-m/sqrt(kappa) integral a^3 A_centered.DS_eta, at fixed g and same all-order CD state",
        "clock_normalized_variance_upper": normalized,
        "canonical_test_force_variance_upper": canonical,
        "unsuppressed_canonical_covariance_coefficient": coeff,
        "clock_normalized_standard_deviation_display": "2e-392 delta ||eta_u||H3 at kappa1e800",
        "canonical_standard_deviation_display": "2e8 delta ||eta_Phi||H3; the kappa suppression cancels in this convention",
        "canonical_dictionary": "eta_Phi=sqrt(kappa) eta_u and Phi_pert=sqrt(kappa) psi; F_Phi[eta_Phi]=F_u[eta_Phi/sqrt(kappa)]",
        "same_full_source_canonical_mean_bound_coefficient": mean,
        "explicit_smaller_ball": SMALL_DELTA,
        "smaller_ball_canonical_source_noise_standard_deviation_upper": 2
        * s.Rational(1, 10**6),
        "smaller_ball_canonical_source_mean_relative_upper": 4 * s.Rational(1, 10**18),
        "small_ball_boundary": "These are source-sector smeared H3 force bounds, not pointwise noise, stress-tensor covariance, a light/metric solution or a complete interacting quantum correction.",
        "checks": checks,
        "gates": {
            "full_clock_normalized_variance_below_display": coeff / K
            < 4 * s.Rational(1, 10**784),
            "canonical_coefficient_below_four_e_sixteen": coeff < 4 * 10**16,
            "smaller_ball_in_full_declared_clock_domain": 0
            < SMALL_DELTA
            < derivatives.DELTA,
            "smaller_ball_canonical_noise_bound": coeff * SMALL_DELTA**2
            < (2 * s.Rational(1, 10**6)) ** 2,
            "same_prior_full_source_mean_estimate": previous.data()[
                "normalized_full_causal_light_force_constant"
            ]
            < 4 * 10**10,
        },
    }
