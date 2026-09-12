"""Uniform Fourier-Laplace bound and ordered bounded-channel inversion."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes

from . import analytic

M = s.Symbol("bounded_channel_norm", nonnegative=True)
SIGMA = 2 * modes.MASS * s.exp(32 + 8 * M)
GAP = (32 + 13 * M) / 5
KBOUND = 5 / (32 + 13 * M)
INVERSE_BOUND = 5 / (32 + 8 * M)


def size(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Rational)):
        raise TypeError("Require an exact finite nonnegative channel norm bound")
    value = s.Rational(value)
    if value < 0:
        raise ValueError("Require a nonnegative channel norm bound")
    return value


def profile(value):
    value = size(value)
    return {
        "channel_norm_bound": value,
        "Laplace_weight": SIGMA.subs(M, value),
        "factor_gap": GAP.subs(M, value),
        "original_channel_inverse_bound": KBOUND.subs(M, value),
        "corrected_channel_inverse_bound": INVERSE_BOUND.subs(M, value),
    }


@cache
def data():
    sigma = s.Symbol("positive_Laplace_weight", positive=True)
    omega = s.Symbol("time_Fourier_frequency", real=True)
    q = s.Symbol("squared_spatial_transfer", nonnegative=True)
    p = (sigma + s.I * omega) ** 2 + q
    modulus = s.expand(s.re(p) ** 2 + s.im(p) ** 2)
    defect = (q - omega**2) ** 2 + sigma**2 * (q + omega**2)
    scalar_identity = s.factor(1 / (GAP - M) - INVERSE_BOUND)
    checks = {
        "full_complex_line_joint_frequency_identity": s.expand(
            modulus - sigma**2 * (sigma**2 + omega**2 + q) - defect
        ),
        "first_sheet_no_cut_on_positive_weight_line": s.im(p) - 2 * sigma * omega,
        "real_frequency_zero_gives_positive_radial_point": p.subs(omega, 0)
        - sigma**2
        - q,
        "chosen_weight_log_radius": s.simplify(
            s.log(SIGMA**2 / (4 * modes.MASS**2)) - 64 - 16 * M
        ),
        "actual_common_gap_at_chosen_weight": s.factor(
            analytic.COEFFICIENT * (64 + 16 * M) - analytic.OFFSET - GAP
        ),
        "original_diagonal_inverse_bound": s.factor(KBOUND * GAP - 1),
        "uniform_contraction_strict_margin": s.factor(
            s.Rational(5, 13) - M * KBOUND - 160 / (13 * (32 + 13 * M))
        ),
        "corrected_gap_retains_full_bounded_term": s.factor(GAP - M - (32 + 8 * M) / 5),
        "ordered_Neumann_sum_bound": scalar_identity,
        "Neumann_bound_before_simplifying": s.factor(
            KBOUND / (1 - M * KBOUND) - INVERSE_BOUND
        ),
    }
    return {
        "joint_frequency_lower_bound": "For lambda=sigma+i omega, q>=0: |lambda^2+q|^2>=sigma^2(sigma^2+omega^2+q), hence |p|>=sigma^2. The bound is uniform over both time frequency and spatial momentum.",
        "selected_positive_weight": SIGMA,
        "uniform_common_real_part_gap": GAP,
        "original_weighted_channel_inverse_bound": KBOUND,
        "bounded_correction_contraction_bound": M * KBOUND,
        "corrected_channel_inverse_bound": INVERSE_BOUND,
        "spaces": "H_sigma^r=L2([0,infinity),exp(-2sigma t)dt;H_x^r(C^2)), every realr. The original causal Kdiag has multiplier diag(1/Ftrace,(3/8)/F2). Plancherel on each Fourier-Laplace line gives its uniform bound. Causal extension preserves support; no pole is removed.",
        "admissible_correction": "V(t) is a measurable time-dependent2x2 channel matrix with essential-sup Euclidean norm<=M. It may mix trace and shear and need not commute with Kdiag. Extend it by zero beyond the requested finite slab. This is not an assertion about the actual unmatched curved response.",
        "ordered_inverse": "K_V=(I+Kdiag V)^-1 Kdiag=Kdiag(I+V Kdiag)^-1. Both converge causally since ||Kdiag V||,||V Kdiag||<=5M/(32+13M)<5/13<1/2. Both inverse identities for Fdiag+V retain the complete distributional initial boundary.",
        "graph": "u in H_sigma^r with (Fdiag+V)u an ordinary H_sigma^r source as a causal distribution, including initial atoms. This graph has the specified unique bounded inverse, not a bounded forward self-map.",
        "checks": checks,
        "gates": {
            "positive_common_gap": GAP.is_positive,
            "strict_positive_corrected_gap": (GAP - M).is_positive,
            "universal_contraction_ratio_below_half": s.Rational(5, 13)
            < s.Rational(1, 2),
            "mass_fixed_at_original_value": modes.MASS == 1000,
            "chosen_weight_above_both_curvature_exponents": 2 * modes.MASS > 108,
            "shear_diagonal_scaling_does_not_worsen_bound": s.Rational(3, 8) < 1,
        },
    }
