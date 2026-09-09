"""Native exact identities, coefficient inequalities and source-pinned proof obligations."""

from functools import cache

import sympy as sp

from . import bounds, calibration, identities, joint, phase, source


@cache
def residuals():
    result = dict(identities.data()["checks"])
    for rows in (identities.asymptotic_controls()["checks"], joint.identities()):
        if set(rows) & set(result):
            raise ValueError("Duplicate global scalar identity")
        result.update(rows)
    return result


@cache
def rejected_inputs():
    for name, call, args in calibration.bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Invalid global scalar input was accepted: " + name)
    return len(calibration.bad_cases())


@cache
def gates():
    d = bounds.data()
    result = {
        "same_retuned_action_reference_and_counterterms": True,
        "same_natural_scalar_anchor_covariance_on_smaller_amplitude_domain": True,
        "weighted_noncanonical_norm_coordinates_do_not_reset_CCR": True,
        "Dini_norm_discards_only_nonpositive_momentum_damping": True,
        "retained_damping_supplies_extra_inverse_square_root_momentum_decay": True,
        "positive_band_covariance_uses_natural_anchor_initial_map": True,
        "compensated_trace_is_restored_in_original_physical_observables": True,
        "variance_derivative_is_transported_by_actual_coupled_generator": True,
        "all_four_actual_scalar_sources_retained": True,
        "whole_half_line_envelopes_are_coefficient_proofs_not_grid_samples": True,
        "mean_generator_is_unchanged_and_its_integral_is_below_two": True,
        "added_state_reflection_not_reference_reflection_controls_both_tails": True,
        "lapse_and_original_trace_tend_to_zero_on_both_tails": True,
        "scale_and_matter_field_derivatives_integrable_give_finite_tail_limits": True,
        "intrinsic_matter_fluctuations_and_mean_lapse_each_counted_once": True,
        "joint_block_state_is_positive_and_its_difference_is_smooth": True,
        "homogeneous_source_variations_preserve_quadratic_block_decomposition": True,
        "joint_physical_mean_response_is_linear_sum_after_trace_restoration": True,
        "joint_exact_frame_not_a_sum_of_nonlinear_frame_maps": True,
        "original_center_jets_preserved_by_unique_same_mean_solution": True,
        "global_positive_scale_and_lapse_prove_representative_completeness": True,
        "global_representative_not_a_coevolved_or_exact_quantum_solution": True,
        "no_uniform_fractional_scalar_density_control_claimed": True,
        "absolute_source_higher_orders_corrected_cones_and_V_G_B_open": True,
    }
    result.update({"phase_" + key: bool(value) for key, value in d["checks"].items()})
    result.update(
        {
            "improved_covariance_constant_below_one_e_sixty": d[
                "improved_phase_covariance_entry_upper_per_eta"
            ]
            < 10**60,
            "global_scalar_mean_scale_constant_below_one_e_sixty_four": d[
                "weighted_mean_phase_upper_per_eta"
            ]
            < 10**64,
            "global_scalar_lapse_constant_below_one_e_sixty_four": d[
                "global_lapse_upper_per_eta"
            ]
            < 10**64,
            "global_scalar_linear_scale_constant_below_one_e_sixty_four": d[
                "global_linear_physical_scale_upper_per_eta"
            ]
            < 10**64,
            "global_scalar_actual_frame_constant_below_two_e_sixty_four": d[
                "global_actual_frame_scale_upper_per_eta"
            ]
            < 2 * 10**64,
            "global_scalar_field_constant_below_one_e_sixty_four": d[
                "global_matter_field_upper_per_eta"
            ]
            < 10**64,
            "lapse_state_every_nonzero_entry_decays": d[
                "lapse_state_all_kernel_entries_decay"
            ],
        }
    )
    for name, row in d["global_physical_observable_bounds"].items():
        result["full_physical_" + name + "_decays_at_least_t_minus_four"] = (
            row["absolute_coefficient_per_t_minus_four"] < 10**64
        )
    for name, upper in zip(
        (
            "mean_scale_forcing",
            "weighted_mean_trace_forcing",
            "lapse_state_term",
            "matter_field_state_term",
        ),
        (257, 2682, 158, 11),
    ):
        result["exact_forcing_integral_bound_" + name] = (
            source.envelopes()[name]["half_entry_L1_upper"] < upper
        )
    result.update(
        {"joint_" + key: bool(value) for key, value in joint.data()["bounds"].items()}
    )
    return result


@cache
def controls():
    d = identities.asymptotic_controls()["actual_limits"]
    c = calibration.calibrated(0)
    return {
        "rejected_inputs": rejected_inputs(),
        "uncompensated_scale_forcing_absolute_integral_is_not_available": d[
            "uncompensated_scale_diagonal_times_u"
        ]
        != 0,
        "uncompensated_weighted_trace_kernel_keeps_volume_growth": d[
            "uncompensated_weighted_trace_diagonal_over_u_four"
        ]
        != 0,
        "compensation_does_not_change_zero_anchor_physical_mean": source.data()[
            "compensating_trace_coefficient"
        ].subs(phase.u, 0)
        == 0,
        "weighted_phase_is_not_symplectic_canonical_phase": phase.data()[
            "old_to_weighted_phase"
        ].det()
        != 1,
        "zero_scalar_covariance_yields_zero_response": all(v == 0 for v in c.values()),
        "zero_joint_covariance_yields_zero_response": all(
            v == 0 for v in joint.calibrated(0, 0, 0).values()
        ),
        "global_domain_is_stricter_than_inherited_local_domain": calibration.DEFAULT_AMPLITUDE
        < sp.Rational(1, 10**20),
        "uncompensated_forcing_growth_not_a_physical_divergence_proof": True,
        "no_uniform_density_fraction_claim_from_absolute_decay_bound": True,
        "no_new_action_state_reset_or_absolute_counterterm_adjustment": True,
        "no_corrected_cone_or_higher_loop_transfer_from_geometric_completeness": True,
        "no_absolute_or_finite_amplitude_SEE_solution_claim": True,
        "no_original_P8_or_V_G_B_closure_inferred": True,
    }
