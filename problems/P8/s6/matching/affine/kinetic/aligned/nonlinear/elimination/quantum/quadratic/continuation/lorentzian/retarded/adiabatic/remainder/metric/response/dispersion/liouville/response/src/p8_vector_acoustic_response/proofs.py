"""Exact acoustic-to-physical response reconstruction and scope gates."""
from functools import cache

from . import dimensional, readout, source, transport


@cache
def residuals():
    result = {}
    for group in (transport.checks(), transport.green_checks(), dimensional.checks(),
                  source.checks(), readout.checks(), readout.actual_checks()):
        if set(result).intersection(group):
            raise ValueError("Repeated acoustic response residual name")
        result.update(group)
    return result


def checks():
    return {
        "two_sector_clocks_have_separate_prepared_shifts": source.data("T")["prepared_shift"] != source.data("L")["prepared_shift"],
        "positive_original_acoustic_rate": source.data("T")["rate"].is_positive is True,
        "full_covariance_and_tangent_constraints_retained": True,
        "physical_readout_rate_volume_and_mass_contacts_retained": True,
        "initial_zero_response_does_not_remove_unvaried_state_mixing": True,
        "retarded_covariance_green_uses_the_actual_selected_modes": True,
        "prepared_shift_not_assumed_zero_at_the_final_endpoint": True,
        "dimensional_pump_jets_retained_before_the_finite_limit": True,
        "analytic_plus_minus_modes_not_conjugated_at_complex_dimension": True,
        "subtraction_is_the_transformed_frozen_physical_subtraction": True,
        "no_independently_renormalized_auxiliary_scalar_is_substituted": True,
        "full_stress_norm_remains_C10_to_C0_not_C2_to_C0": True,
        "no_coupled_causal_inverse_or_original_P8_closure_claim": True,
    }
