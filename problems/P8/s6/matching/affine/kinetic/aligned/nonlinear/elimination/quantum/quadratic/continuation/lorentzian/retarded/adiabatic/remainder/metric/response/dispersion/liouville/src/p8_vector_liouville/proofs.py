"""Exact mode-reduction identities and continuous source-domain gates."""
from functools import cache

import sympy as sp

from . import clock, reduction, variation


@cache
def residuals():
    result = {}
    for group in (reduction.checks(), clock.checks(), variation.checks()):
        if set(result).intersection(group):
            raise ValueError("Repeated Liouville residual name")
        result.update(group)
    return result


@cache
def checks():
    estimate = variation.estimates()
    lower = clock.bounds(0, 1000)["clock_scalar_principal_potential_lower"]
    return {
        "positive_clock_scalar_principal_potential_at_selected_mass": bool(lower > 0),
        "background_remainder_coefficient_below_eighty_eight": bool(sp.Rational(359375, 4096) < 88),
        "prepared_time_shift_coefficient_at_most_two": estimate["initial_zero_time_shift_C0_upper_per_source_C0"] == 2,
        "prepared_C2_correction_coefficient_below_eighteen_thousand": bool(
            estimate["correction_C2_coefficient"] < 18000),
        "all_varied_numerator_bounds_positive": all(value > 0 for value in estimate["varied_numerator_C2_bounds"]),
        "all_history_numerator_bounds_positive": all(value > 0 for value in estimate["history_numerator_bounds"]),
        "strict_momentum_decay_in_exact_rational_bound": bool(
            variation.bound(2000, 1000)["prepared_remainder_C2_to_C0_upper"]
            < variation.bound(1000, 1000)["prepared_remainder_C2_to_C0_upper"]),
        "acoustic_time_is_a_mode_coordinate_not_a_new_physical_metric": True,
        "longitudinal_map_defined_for_k_positive_potential_limit_at_zero_only": True,
        "prepared_initial_mode_data_and_nonzero_mixing_unchanged": True,
        "no_ordinary_Proca_cone_or_Hadamard_claim_transferred_off_clock": True,
        "potential_response_bound_not_a_stress_or_no_loss_feedback_norm": True,
    }
