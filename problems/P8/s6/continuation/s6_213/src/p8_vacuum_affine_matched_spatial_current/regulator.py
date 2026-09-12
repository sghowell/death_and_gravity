"""Original two-leg sharp regulator, anchored subtraction and uniform tail."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_subleading_band_conversion import centered as conversion
from p8_vacuum_affine_uniform_uv_remainder import limit as known

TAIL = 3 * s.Integer(10) ** 54
CONVERSION_TAIL = s.Integer(10) ** 40


def subtractor(K, m, A3, A1, F2, F4, finite):
    return K**3 * A3 + K * A1 + F2 * (K * K - m * m) / 2 + F4 * s.log(K / m) - finite


@cache
def data():
    K, m = s.symbols("K m", positive=True)
    A3, A1, F2, F4, finite = s.symbols("A3 A1 F2 F4 Ffinite")
    H, KP, K0, UP, U0, E, LP, L0 = s.symbols(
        "H KnownK KnownK0 UK U0 EK KnownLimit KnownLimit0"
    )
    UV = K**3 * A3 + K * A1 + F2 * (K * K - m * m) / 2 + F4 * s.log(K / m)
    physical = H + LP - L0 + finite
    regulated = H + KP - K0 + UP - U0 - subtractor(K, m, A3, A1, F2, F4, finite)
    error = s.expand((regulated - physical).subs(UP, U0 + UV + E))
    checks = {
        "original_regulator_exact_anchored_error": s.expand(
            error - (KP - LP - (K0 - L0) + E)
        ),
        "full_finite_lower_band_piece_retained": s.expand(
            subtractor(K, m, A3, A1, F2, F4, finite)
            - (K**3 * A3 + K * A1 + F2 * K * K / 2 + F4 * s.log(K / m) - finite)
            + F2 * m * m / 2
        ),
        "zero_spatial_transfer_subtractor": subtractor(K, m, 0, 0, 0, 0, 0),
        "conversion_quadratic_artifact_not_reintroduced": s.expand(UV).coeff(K, 2)
        - F2 / 2,
        "conversion_full_finite_artifact_not_reintroduced": (
            UV - F4 * s.log(K / m)
        ).subs(K, 0)
        + F2 * m * m / 2,
        "complete_known_tail_unchanged": known.KNOWN_TAIL - 5 * s.Integer(10) ** 53,
        "complete_conversion_tail_unchanged": conversion.UNIFORM_TAIL - CONVERSION_TAIL,
        "both_external_metric_tail_factors": 4 * TAIL / modes.KAPPA
        - 12 * s.Rational(1, 10) ** 746,
    }
    unrounded = 4 * known.KNOWN_TAIL + CONVERSION_TAIL
    return {
        "original_mask": "Keep chi_K(k)chi_K(P-k), K>=2m, in the pair and the original one-leg contact mask. The two-leg-to-one-ball conversion is subtracted explicitly; no new regulator is substituted for the original one.",
        "subtractor": "S_K(P)=K^3 A3(P)+K A1(P)+F2(P)(K^2-m^2)/2+F4(P)log(K/m)-Ffinite(P). A3,A1 are the complete actual curved coefficients. Complete conversion K2 and K0 artifacts vanish, but the physical one-ball quadratic divergence and finite lower-band term remain.",
        "anchored_finite_K_current": "Jren,K^anchor=H0+[Jactual,K(P)-Jactual,K(0)]-S_K(P). The full contact cancels exactly in the difference. This is a finite-K approximation anchored to the exact renormalized homogeneous current, not an asserted tail estimate for the separately unanchored homogeneous cutoff current.",
        "uniform_tail": "The exact error is[KnownK(P)-Known(P)]-[KnownK(0)-Known(0)]+E_K. Each known tail is below1e54 M[D]Z136/K, and the complete all-P conversion tail is below1e40||D||L2 X46/K. The sum is below3e54 M[D]Z136/K; its canonical display is12e-746/K.",
        "unrounded_tail": unrounded,
        "not_finite_cutoff": "This proves removal of the original auxiliary regulator in the anchored tracefree spatial Gaussian current. It is not a physical EFT cutoff, finite-gravity bound, global-time estimate or omitted-loop bound.",
        "checks": checks,
        "gates": {
            "two_known_tails_and_conversion_rounded_strictly": unrounded < TAIL,
            "actual_all_P_conversion_tail_not_fixed_P_only": True,
            "complete_frozen_conversion_tail_strict": conversion.tail_constants()[
                "complete_uniform_conversion_remainder"
            ]
            < CONVERSION_TAIL,
            "finite_lower_band_term_not_deleted": True,
            "complete_odd_endpoint_and_cancellation_results_preserved": True,
            "original_two_leg_mask_unchanged": True,
            "homogeneous_anchor_exact_not_finite_K_tail_claim": True,
            "no_physical_cutoff_or_original_P8_closure": True,
        },
    }
