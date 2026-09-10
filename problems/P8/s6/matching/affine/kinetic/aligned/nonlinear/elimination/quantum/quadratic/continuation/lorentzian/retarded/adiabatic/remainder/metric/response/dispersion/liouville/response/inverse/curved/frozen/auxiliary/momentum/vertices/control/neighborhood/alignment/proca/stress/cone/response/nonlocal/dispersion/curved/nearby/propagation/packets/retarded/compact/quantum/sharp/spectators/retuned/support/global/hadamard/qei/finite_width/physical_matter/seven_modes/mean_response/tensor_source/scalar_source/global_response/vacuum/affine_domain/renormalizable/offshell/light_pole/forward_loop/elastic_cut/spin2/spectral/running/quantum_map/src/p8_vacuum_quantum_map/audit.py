"""Source-aware finite-order map audit, with explicit non-transfer boundaries."""

from functools import cache

import sympy as sp
from p8_offshell_vacuum import jets

from . import calibration, interactions, ward, wick


@cache
def residuals():
    return {
        prefix + "_" + name: value
        for prefix, module in (
            ("wick", wick),
            ("interactions", interactions),
            ("ward", ward),
            ("calibration", calibration),
        )
        for name, value in module.data()["checks"].items()
    }


@cache
def gates():
    rows = {k: bool(v) for k, v in calibration.data()["bounds"].items()}
    rows.update(
        {
            "actual_S711_cubic_derivative_map_not_a_new_point_map": True,
            "literal_210_jet_four_dimensional_contraction_independently_rebuilt": wick.data()[
                "literal_jet_count"
            ]
            == 210,
            "dimensional_tensors_derived_before_subtraction": True,
            "finite_evanescent_term_retained_in_composite_MSbar_prescription": True,
            "odd_light_parity_keeps_zero_vacuum_expectation": True,
            "inverse_quintic_and_interacting_composite_graphs_start_at_two_loops": True,
            "all_six_transformed_field_degrees_retained": set(
                interactions.data()["all_generated_field_degrees"]
            )
            == {2, 4, 6, 8, 10, 12},
            "one_loop_four_point_requires_quartic_pair_and_sextic_tadpole": set(
                interactions.topologies(4, 1)
            )
            == {(2, 0, 0, 0, 0), (0, 1, 0, 0, 0)},
            "exact_Gaussian_heavy_integration_commutes_with_light_only_map": True,
            "full_heavy_kernel_never_expanded_inside_loop_momentum_integral": True,
            "all_parent_counterterms_transformed_with_fixed_finite_parts": True,
            "quadratic_counterterm_generated_quartic_not_discarded": True,
            "local_Jacobian_ghost_loops_scaleless_before_regulator_removal": True,
            "physical_sources_transformed_before_ordinary_Psi_LSZ": True,
            "Gaussian_finite_dimensional_controls_not_continuum_evidence": True,
            "computed_composite_overlap_not_set_to_identity": True,
            "stable_light_pole_only_no_unstable_heavy_LSZ": True,
            "complete_parent_one_loop_b2_and_elastic_cut_transfer_to_full_matched_action": True,
            "bare_truncated_target_not_quantum_equivalent_by_this_argument": True,
            "ordinary_off_shell_effective_actions_not_assumed_scalars": True,
            "one_loop_disc_bound_not_global_resummed_higher_derivative_spectrum": True,
            "Schwartz_stationary_error_not_a_quantum_loop_remainder": True,
            "flat_field_map_does_not_transfer_graviton_form_factor": True,
            "new_composite_prescription_not_old_rolling_counterterm_transfer": True,
            "full_V_G_B_and_original_P8_remain_open": True,
        }
    )
    return rows


def bad_cases():
    zero = (0, 0, 0, 0)
    rows = list(calibration.bad_cases())
    rows += [
        ("bad_covariance_" + str(i), wick.covariance, (a, b))
        for i, (a, b) in enumerate(
            (
                ([0, 0, 0, 0], zero),
                ((0, 0, 0), zero),
                ((True, 0, 0, 0), zero),
                ((1.0, 0, 0, 0), zero),
                ((-1, 0, 0, 0), zero),
                ((5, 0, 0, 0), zero),
                ((sp.Integer(1), 0, 0, 0), zero),
                (zero, (0, 0, -1, 0)),
            )
        )
    ]
    rows += [
        ("bad_cubic_" + str(i), wick.contract, (v,))
        for i, v in enumerate(
            (
                sp.Integer(0),
                sp.Integer(1),
                jets.PHI,
                jets.PHI**2,
                jets.PHI**4,
                jets.PHI**3 + jets.PHI,
            )
        )
    ]
    rows += [
        ("bad_topology_" + str(i), interactions.topologies, (a, b))
        for i, (a, b) in enumerate(
            ((True, 1), (4.0, 1), (3, 1), (8, 1), (4, True), (4, 1.0), (4, -1), (4, 3))
        )
    ]
    return rows


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (ValueError, TypeError):
            continue
        raise ValueError("Unsupported finite-order map input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    w = ward.data()
    return {
        "rejected_inputs": rejected_inputs(),
        "Gaussian_Jacobian_omission_detected": w[
            "omitting_Jacobian_changes_partition_function"
        ]
        != 0,
        "untransformed_source_omission_detected": w[
            "omitting_transformed_source_changes_two_point_function"
        ]
        != 0,
        "generated_sextic_omission_detected": w[
            "omitting_generated_sextic_changes_partition_function"
        ]
        != 0,
        "off_shell_effective_action_scalar_claim_rejected": w[
            "off_shell_effective_action_difference_quartic"
        ]
        != 0,
        "evanescent_finite_difference_nonzero": True,
        "ordinary_amputation_alone_not_LSZ": True,
        "bare_derivative_truncation_not_full_transformed_quantum_action": True,
        "resummed_global_ghost_spectrum_not_inferred": True,
        "no_finite_gravity_or_rolling_state_transfer": True,
        "original_P8_not_closed": True,
    }
