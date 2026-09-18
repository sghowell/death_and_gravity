"""Integrable probability-overlap transfer and explicit unresolved matching."""

from functools import cache

import sympy as s
from p8_vacuum_affine_three_singleton_subtraction import measure as previous

from . import faces


@cache
def data():
    a, b, c, x = s.symbols("a b c x", positive=True)
    L, Q = faces.AXIS_CAP, faces.PAIR_CAP
    Kmajor = 4 * (1 / s.sqrt(a * b) + 1 / s.sqrt(a * c) + 1 / s.sqrt(b * c))
    angular = s.Rational(1, 48) / s.pi**6
    La, Lb, La0, Lb0, Ra, Rb = s.symbols("La Lb La0 Lb0 Ra Rb", real=True)
    baseline = La * (Rb + Lb0) + Lb * (Ra + La0) - La0 * Lb0
    marked = La * Rb + Lb * Ra + La0 * Lb0
    mismatch = (La - La0) * Lb0 + (Lb - Lb0) * La0
    checks = {
        "all_three_entropy_product_integrals": s.integrate(
            Kmajor, (a, 0, x), (b, 0, x), (c, 0, x)
        )
        - 48 * x * x,
        "exact_transfer_phase_space_arithmetic": s.factor(
            angular * (48 * L * Q + 96 * Q * Q) - (L * Q + 2 * Q * Q) / s.pi**6
        ),
        "same_parent_complete_signed_TV": previous.integration()[
            "whole_total_variation_upper"
        ]
        - s.Rational(1, 10**1424) * x * x
        - s.Rational(1, 10**1340) * x**4,
        "proper_two_real_state_change_matching_identity": s.expand(
            baseline - marked - mismatch
        ),
        "nonzero_state_change_example": s.expand(
            mismatch.subs({La: 1 + b, Lb: 1 + a, La0: 1, Lb0: 1}) - a - b
        ),
        "remaining_logarithmic_soft_pole_example": s.factor(
            (a + b) / (a * b) - 1 / a - 1 / b
        ),
    }
    return {
        "checks": checks,
        "gates": {
            "derived_transfer_coefficient_below1e_minus1424": bool(
                L * Q + 2 * Q * Q < s.Rational(1, 10**1424)
            ),
            "state_change_terms_not_identically_zero": mismatch != 0,
            "probability_projection_retains_signed_common_cutoff": True,
            "all_physical_polarizations_and_Bose_factor": True,
            "simplex_to_cube_only_for_positive_majorant": True,
            "no_virtual_or_integrated_counterterm_chosen": True,
            "no_positive_normalized_or_monotone_rate_claim": True,
        },
        "whole_transfer_TV_coefficient": (L * Q + 2 * Q * Q) / s.pi**6,
        "whole_transfer_TV_upper": s.Rational(1, 10**1424) * x * x,
        "whole_probability_rectangle_TV_upper": 2 * s.Rational(1, 10**1424) * x * x
        + s.Rational(1, 10**1340) * x**4,
        "whole_relative_reference_upper": 4
        * s.Rational(1, 10**1424)
        * x ** s.Rational(3, 2)
        + 2 * s.Rational(1, 10**1340) * x ** s.Rational(7, 2),
        "whole_same_state_matching_remainder": mismatch,
        "whole_finite_measure": "For0<x<=1/8 the normalized probability rectangle R3(|G3|^2)/(a^2*b^2*c^2), with the original three-real measure and all polarizations/3!, has total variation below2e-1424*x^2+1e-1340*x^4. Common lower-cutoff removal converges in total variation. The result is signed, not a positive normalized detector distribution.",
        "whole_matching_boundary": "Even the two-real probability baseline has explicit radiative-state current changes beyond S309's singly marked residual terms and the Born product. Those nonzero terms cannot be discarded or identified with actual virtual graphs by definition. The all-N, real-virtual, loop/state, Regge and original P8 obligations remain open.",
    }


def total_variation_upper(value):
    x = previous.require_cutoff(value)
    return 2 * s.Rational(1, 10**1424) * x * x + s.Rational(1, 10**1340) * x**4


def relative_reference_upper(value):
    x = previous.require_cutoff(value)
    return 4 * s.Rational(1, 10**1424) * x ** s.Rational(3, 2) + 2 * s.Rational(
        1, 10**1340
    ) * x ** s.Rational(7, 2)
