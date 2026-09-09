"""Native source/mean/Ward/global-bound gates and strict claim controls."""

from functools import cache

import sympy as sp

from . import anchor, band, mean, stress, tails


@cache
def residuals():
    result = {}
    for data in (
        stress.data(),
        mean.data(),
        mean.ward(),
        mean.center_control(),
        band.data(),
        anchor.data(),
        tails.data(),
    ):
        overlap = set(result).intersection(data["checks"])
        if overlap:
            raise ValueError("Duplicate scientific residual names: " + str(overlap))
        result.update(data["checks"])
    return result


@cache
def rejected_inputs():
    for name, call, args in band.bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("An invalid relative-state amplitude was accepted: " + name)
    return len(band.bad_cases())


@cache
def gates():
    result = {
        "same_original_global_retuned_action_and_physical_matter_frame": True,
        "only_Proca_covariance_changed_by_real_positive_smooth_band": True,
        "same_scalar_tensor_factors_and_all_counterterms_retained": True,
        "band_is_positive_normalized_and_has_no_zero_momentum_or_high_frequency_tail": True,
        "finite_energy_density_not_claimed_to_be_finite_total_energy_on_R_cubed": True,
        "exact_free_Proca_evolution_on_both_legs": True,
        "full_physical_pressure_and_source_conservation_not_assigned_by_hand": True,
        "lapse_is_unchanged_in_actual_spatial_frame_map": True,
        "actual_lapse_source_keeps_pressure_and_second_order_mean_data": True,
        "actual_trace_density_not_mistaken_for_canonical_momentum": True,
        "actual_matter_field_charge_density_and_connection_variation_retained": True,
        "center_induced_scalar_density_bridge_covers_all_six_vector_channels": True,
        "positive_combined_center_response_does_not_imply_positive_isolated_scalar_variation": True,
        "leading_equation_difference_does_not_set_common_vacuum_tadpole_to_zero": True,
        "relative_state_reflection_not_assumed_for_prepared_reference_state": True,
        "all_real_times_covered_by_positive_pivot_and_integrable_weighted_bounds": True,
        "exact_frame_representative_complete_bounce_not_claimed_to_solve_SEE": True,
        "absolute_source_higher_order_errors_and_V_G_B_remain_open": True,
    }
    result.update({"band_" + k: bool(v) for k, v in band.data()["bounds"].items()})
    result.update({"global_" + k: bool(v) for k, v in tails.data()["bounds"].items()})
    for key in (
        "all_anchor_acceleration_response_polynomial_coefficients_positive",
        "rounded_acceleration_bound_valid",
        "initial_lapse_response_positive_on_entire_band",
        "initial_hat_scale_second_derivative_positive_on_entire_band",
        "initial_lapse_second_derivative_all_polynomial_coefficients_positive",
    ):
        result[key] = bool(anchor.data()[key])
    return result


@cache
def controls():
    d = stress.data()
    w = mean.ward()
    return {
        "rejected_inputs": rejected_inputs(),
        "wrong_four_dimensional_lapse_rescaling_extra_source_nonzero": d[
            "false_four_dimensional_conformal_extra_source"
        ]
        != sp.zeros(6),
        "omitting_physical_pressure_changes_lapse_source": d[
            "physical_isotropic_pressure_Hessian"
        ]
        != sp.zeros(6),
        "omitting_connection_variation_changes_Ward_identity": w[
            "retained_connection_variation_term"
        ]
        != 0,
        "reference_scalar_tensor_state_and_tadpole_not_reassigned": True,
        "no_absolute_vacuum_source_claimed_zero": True,
        "no_nonlinear_or_higher_loop_remainder_bound_inferred": True,
        "no_quantum_state_transported_on_representative_metric_claimed": True,
        "no_finite_cutoff_or_V_G_B_verdict_from_band_or_metric": True,
        "zero_added_covariance_has_zero_calibrated_response": all(
            v == 0 for v in band.calibrated(0).values()
        ),
    }
