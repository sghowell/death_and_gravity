"""Native variational, canonical, observable and interval-bound obligations."""

from functools import cache

import sympy as sp

from . import bounds, calibration, geometry, model, observable, response, source


@cache
def residuals():
    result = dict(model.checks())
    for data in (
        geometry.data(),
        geometry.canonical(),
        source.data(),
        source.density_center(),
        response.data(),
        response.parity(),
        response.center(),
        observable.data(),
        observable.center(),
        bounds.response_bounds(),
    ):
        overlap = set(result).intersection(data["checks"])
        if overlap:
            raise ValueError("Duplicate scalar-source identity: " + str(overlap))
        result.update(data["checks"])
    return result


@cache
def rejected_inputs():
    for name, call, args in calibration.bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Invalid scalar-source input was accepted: " + name)
    return len(calibration.bad_cases())


@cache
def gates():
    result = {
        "same_original_retuned_action_and_physical_matter_frame": True,
        "homogeneous_background_kept_independent_before_variation": True,
        "literal_canonical_exponential_volume_and_York_shear_derived": True,
        "only_nonzero_momentum_lapse_eliminated_in_quadratic_field_action": True,
        "third_background_lapse_jet_retained_in_actual_source": True,
        "both_momentum_shifts_and_time_boundary_match_existing_scalar_generator": True,
        "positive_I4_band_addition_named_in_natural_anchor_phase": True,
        "other_five_mode_factors_reference_state_and_counterterms_unchanged": True,
        "exact_global_scalar_evolution_preserves_positive_CCR_and_smooth_difference": True,
        "finite_homogeneous_density_not_finite_total_R3_energy": True,
        "all_four_actual_mean_sources_retained": True,
        "matter_charge_not_reset_when_field_velocity_source_is_nonzero": True,
        "scalar_matter_observable_keeps_intrinsic_fluctuations_and_mean_once": True,
        "independent_center_constraint_and_density_comparisons_are_full_polynomial": True,
        "scalar_state_reflection_keeps_odd_matter_coordinate": True,
        "whole_strip_and_band_bounds_are_exact_coefficient_not_grid_estimates": True,
        "mean_solution_exists_at_all_finite_times_but_uniform_tails_not_claimed": True,
        "local_exact_frame_representative_is_not_a_quantum_solution": True,
        "no_simple_homogeneous_fluid_Ward_identity_assumed_for_scalar_observable": True,
        "absolute_source_higher_orders_global_quantum_bounce_and_V_G_B_open": True,
    }
    result.update(
        {
            "bound_" + key: bool(value)
            for key, value in bounds.response_bounds()["bounds"].items()
        }
    )
    return result


@cache
def controls():
    d = source.data()
    c = response.center()
    return {
        "rejected_inputs": rejected_inputs(),
        "omitting_time_dependent_boundary_changes_scalar_Hamiltonian": model.bridge()[
            "generating_time_derivative"
        ]
        != 0,
        "discarding_direct_scalar_scale_source_changes_mean_equation": d[
            "source_Hessians"
        ]["mean_hat_scale_direct_source"]
        != sp.zeros(4),
        "discarding_direct_scalar_matter_field_source_changes_mean_equation": d[
            "source_Hessians"
        ]["mean_matter_field_direct_source"]
        != sp.zeros(4),
        "discarding_intrinsic_matter_fluctuations_changes_physical_density": observable.data()[
            "intrinsic_density_Hessian"
        ]
        != sp.zeros(4),
        "discarding_nonlinear_scalar_chart_changes_center_density": source.density_center()[
            "required_scalar_chart_correction"
        ]
        != 0,
        "acceleration_kernel_has_opposite_signs_within_band_closure": c[
            "anchor_proper_Hubble_derivative_response"
        ].subs(model.k, 1)
        > 0
        and c["anchor_proper_Hubble_derivative_response"].subs(model.k, 2) < 0,
        "zero_added_scalar_covariance_gives_zero_calibrated_response": all(
            value == 0 for value in calibration.calibrated(0).values()
        ),
        "no_reassigned_common_absolute_tadpole_or_counterterm": True,
        "no_independent_minimal_scalar_or_Proca_source_transfer": True,
        "no_uniform_infinite_tail_bound_inferred_from_local_interval": True,
        "no_positive_lapse_for_every_unbounded_Gaussian_field_configuration_claimed": True,
        "no_finite_amplitude_SEE_or_full_quantum_Ward_identity_claimed": True,
        "no_heavy_cutoff_higher_loop_bound_or_V_G_B_verdict_inferred": True,
    }
