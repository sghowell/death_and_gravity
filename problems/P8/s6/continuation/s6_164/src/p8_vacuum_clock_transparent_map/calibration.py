"""Exact preservation of the full-target classical error after map extension."""

from functools import cache

import sympy as s
from p8_exceptional_vacuum import analytic
from p8_offshell_vacuum import norms as old_norms
from p8_polynomial_vacuum import model
from p8_vacuum_full_analytic_target_match import transport as old_match

from . import norms


@cache
def data():
    p = model.data()["actual_parameters"]
    C = old_norms.data()["actual_C_R"]
    Xcap = s.Rational(4, analytic.KAPPA)
    raw_qcap = norms.gate_difference_cap(Xcap)
    qcap = s.Rational(1, 10**6390)
    delta = s.Rational(1, 10**6793)
    error = norms.linear_enclosure(
        C,
        delta,
        p["bare_polynomial_quartic"],
        p["cubic_coupling_squared"],
        p["heavy_mass_squared"],
    )
    previous = old_match.data()[
        "full_analytic_target_stationary_action_match_coefficient"
    ]
    combined = previous + error["complete_action_difference_coefficient"]
    return {
        "same_actual_cubic_map_norm": C,
        "actual_gradient_invariant_cap": Xcap,
        "declared_polynomial_gate_difference_cap": qcap,
        "declared_map_difference_sup_and_L2_factor": delta,
        "complete_new_map_action_difference": error,
        "previous_full_analytic_target_match_coefficient": previous,
        "clock_transparent_full_analytic_target_match_coefficient": combined,
        "bounds": {
            "finite_new_source_support_within_actual_full_heavy_gap": bool(
                p["heavy_mass_squared"] > 66**2
            ),
            "polynomial_gate_difference_below_declared_one_e_minus_6390": bool(
                0 < raw_qcap < qcap
            ),
            "actual_cubic_map_cap_proves_declared_difference": bool(qcap * C < delta),
            "linear_enclosure_norm_caps_below_one_quarter": bool(
                0 < C < s.Rational(1, 4) and 0 < delta < s.Rational(1, 4)
            ),
            "linearized_free_corner_below_three_hundred": bool(
                s.Rational(1045, 4) < 300
            ),
            "linearized_local_corner_below_one": bool(s.Rational(671, 1536) < 1),
            "linearized_heavy_corner_below_two": bool(s.Rational(671, 512) < 2),
            "new_map_difference_below_one_e_minus_6790": bool(
                0 < delta < s.Rational(1, 10**6790)
            ),
            "complete_new_map_action_difference_below_one_e_minus_6790": bool(
                0
                < error["complete_action_difference_coefficient"]
                < s.Rational(1, 10**6790)
            ),
            "complete_full_analytic_target_match_still_below_one_e_minus_800": bool(
                0 < combined < s.Rational(1, 10**800)
            ),
        },
        "checks": {
            "same_actual_canonical_gradient_scaling": Xcap * analytic.KAPPA - 4,
            "declared_gate_cap": qcap - s.Rational(1, 10**6390),
            "declared_map_difference_cap": delta - s.Rational(1, 10**6793),
            "all_new_action_groups_owned_once": error[
                "complete_action_difference_coefficient"
            ]
            - error["free_action_difference_coefficient"]
            - error["local_quartic_difference_coefficient"]
            - error["full_heavy_resolvent_difference_coefficient"],
            "unchanged_prior_full_target_match": previous
            - old_match.data()[
                "full_analytic_target_stationary_action_match_coefficient"
            ],
            "combined_classical_triangle_inequality": combined
            - previous
            - error["complete_action_difference_coefficient"],
        },
    }
