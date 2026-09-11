"""Compose with the previously proved full-resolvent stationary action."""

from functools import cache

import sympy as s
from p8_offshell_vacuum import band
from p8_polynomial_vacuum import calibration as parent_parameters

from . import bounds, calibration, holomorphic


@cache
def data():
    previous = band.data()
    h = holomorphic.data()
    old = previous["combined_tree_action_error_coefficient"]
    added = calibration.data()[
        "full_target_action_minus_quartic_target_common_class_coefficient"
    ]
    combined = bounds.compose(old, added)
    p = parent_parameters.point()
    Xcap = s.Rational(4, h["fixed_kappa"])
    return {
        "classical_common_function_class": previous["domain"],
        "previous_quartic_target_action_error_coefficient": old,
        "full_analytic_target_correction_coefficient": added,
        "full_analytic_target_stationary_action_match_coefficient": combined,
        "new_field_Fourier_radius": previous["new_field_fourier_radius"],
        "mapped_field_Fourier_radius": previous["mapped_field_fourier_radius"],
        "mapped_quadratic_source_Fourier_radius": previous[
            "mapped_quadratic_source_fourier_radius"
        ],
        "centered_box_multiplier_radius": previous["centered_box_spectral_radius"],
        "rank_regular_target_X_absolute_upper_on_real_common_class": Xcap,
        "conclusion": "For the same real Schwartz Psi with Fourier support in the Euclidean four-momentum unit ball and Fourier L1 norm <=1, let U=||Psi||_2. The full polynomial classical action evaluated at Phi=Psi+R[Psi] and the stationary full heavy resolvent differs from the FULL finite-kappa S6.109 analytic flat scalar target action evaluated at Psi by at most the declared combined coefficient times U^2, strictly below 1e-800 U^2. Fermion and gauge classical fields are set to zero. The physical metric is fixed and flat.",
        "bounds": {
            "unchanged_mapped_source_strictly_inside_heavy_resolvent_gap": bool(
                previous["centered_box_spectral_radius"] < p["D"]
            ),
            "real_common_class_inside_S6_109_rank_regular_X_domain": bool(
                Xcap < s.Rational(1, 4 * h["fixed_n"])
            ),
            "full_analytic_target_common_class_match_strictly_below_one_e_minus_800": bool(
                0 < combined < s.Rational(1, 10**800)
            ),
            "target_correction_smaller_than_previous_match_allowance": bool(
                0 < added < old
            ),
        },
        "checks": {
            "unchanged_original_quartic_target_error": old
            - band.data()["combined_tree_action_error_coefficient"],
            "full_target_triangle_inequality": combined - old - added,
            "same_mapped_field_support": previous["mapped_field_fourier_radius"] - 3,
            "same_mapped_heavy_source_support": previous[
                "mapped_quadratic_source_fourier_radius"
            ]
            - 6,
            "same_centered_operator_spectral_radius": previous[
                "centered_box_spectral_radius"
            ]
            - 38,
            "actual_canonical_to_dimensionless_X": Xcap * h["fixed_kappa"] - 4,
        },
        "not_inferred": "Neither a quantum matching theorem for the full analytic target nor a global inverse field map, rolling-state dictionary, cutoff, metric/stress error, finite-gravity or common-parent completion follows. The analytic coefficient circle controls field degree at fixed jets, not integration over unbounded loop momenta.",
    }
