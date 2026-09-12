"""All-momentum endpoint remainder and the original removed two-leg union."""

from functools import cache

import sympy as s
from p8_vacuum_affine_reference_spatial_analyticity import boundary

from . import near, real

FULL = 3 * s.Integer(10) ** 54
TAIL = 5 * s.Integer(10) ** 60


@cache
def constants():
    c = near.constants()
    low = real.constants()["unexpanded_low_band_integral_upper"]
    C4 = c["complete_near_P4_upper"]
    C5 = c["complete_near_P5_upper"]
    B = c["internal_near_radius_over_external_P"] + 1
    return {
        "low_band_upper": low,
        "near_fourth_power_upper": C4,
        "near_fifth_power_upper": C5,
        "removed_near_external_threshold_denominator": B,
        "full_low_near_far_upper": low + C4 + C5 + boundary.FAR,
        "full_removed_union_tail_upper": (real.MASS + 1) * low
        + B * (C4 + C5)
        + boundary.TAIL,
    }


@cache
def data():
    P = s.Symbol("P", nonnegative=True)
    K, d = s.symbols("K d", positive=True)
    overlap = s.pi * (4 * K + d) * (2 * K - d) ** 2 / 12
    c = constants()
    checks = {
        "fourth_power_six_derivative_weight": s.expand(
            (1 + P * P) ** 3 - P**4 - (P**6 + 2 * P**4 + 3 * P**2 + 1)
        ),
        "fifth_power_six_derivative_weight": s.expand(
            (1 + P * P) ** 3
            - P**5
            - (
                P**4 * (P - s.Rational(1, 2)) ** 2
                + s.Rational(11, 4) * P**4
                + 3 * P**2
                + 1
            )
        ),
        "sixth_power_same_weight": s.expand(
            (1 + P * P) ** 3 - P**6 - (3 * P**4 + 3 * P**2 + 1)
        ),
        "first_power_same_weight": s.expand(
            (1 + P * P) ** 3
            - P
            - (
                (P - s.Rational(1, 2)) ** 2
                + s.Rational(3, 4)
                + 2 * P**2
                + 3 * P**4
                + P**6
            )
        ),
        "complete_near_tail_same_spatial_regularity": s.expand(
            P
            * (c["near_fourth_power_upper"] * P**4 + c["near_fifth_power_upper"] * P**5)
            - c["near_fourth_power_upper"] * P**5
            - c["near_fifth_power_upper"] * P**6
        ),
        "constant_integrand_two_ball_intersection": s.expand(
            overlap - (4 * s.pi * K**3 / 3 - s.pi * K * K * d + s.pi * d**3 / 12)
        ),
        "right_derivative_of_moving_band_volume": s.diff(overlap, d).subs(d, 0)
        + s.pi * K * K,
    }
    return {
        "remainder_definition": "At the original two-created-mode regulator, Q_K integrates E(P,k) for |k|<m and E(P,k)-T4E(P,k) for |k|>=m. The high band is partitioned into FAR and NEAR only for estimation. No regulator domain is Taylor-differentiated.",
        "complete_bound": "The unexpanded low part, near C4|P|^4+C5|P|^5 part and S202 far part give |Q_K|<3e54||D||L2 X46[Gamma] uniformly. The same bound holds for its absolute weak continuum limit.",
        "near_removed_union": "On NEAR, |k|<B0|P|, B0=2Amax/delta, and |l|<=(B0+1)|P|. A removed pair therefore forces |P|>K/(B0+1). Multiply the near majorant by(B0+1)|P|/K; the resulting fifth and sixth powers are still covered by X46.",
        "low_removed_union": "On |k|<m, max(|k|,|l|)<m+|P|. The removed-pair indicator is bounded by(m+|P|)/K. The same spatial weight gives the displayed(m+1) low coefficient. This holds for all K>=m, including the original threshold.",
        "complete_tail": "Add the unchanged S202 far removed-union bound. The full tail is below5e60||D||L2 X46[Gamma]/K. Thus Q_K converges absolutely as a weak spatial/time-jet functional, retaining both created momenta.",
        "no_locality_inference": "The subtracted high-band Taylor-integrand sector is not claimed to be polynomial after the moving band intersection is integrated. Its divergent and finite content and matching with the distinct one-leg contact remain outside Q.",
        "geometric_counterexample": "For the constant polynomial integrand1, the intersection of two radiusK balls separated by d=|P| has volume4pi K^3/3-pi K^2|P|+pi|P|^3/12 for |P|<=2K. Its opposite one-sided derivatives at zero differ. This demonstrates why integrated sharp-band polynomiality cannot be assumed; it does not compute an actual Proca divergence coefficient.",
        "constants": c,
        "checks": checks,
        "gates": {
            "complete_low_near_far_display": c["full_low_near_far_upper"] < FULL,
            "complete_regulator_tail_display": c["full_removed_union_tail_upper"]
            < TAIL,
            "same_six_spatial_derivatives_for_tail": max(j + 1 for j in (4, 5)) <= 6,
            "original_both_leg_near_threshold": c[
                "removed_near_external_threshold_denominator"
            ]
            > 1,
            "no_spatial_derivative_of_regulator_domain": True,
            "constant_integrand_sharp_band_derivative_jump": 2 * s.pi * K * K > 0,
        },
    }
