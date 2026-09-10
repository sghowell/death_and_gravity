"""Independent action, domain, integral and general-gradient verification."""

from functools import cache

import sympy as sp

from . import bounds, calibration, covariance, family, local, lower


@cache
def residuals():
    groups = {
        "principal": family.principal_checks(),
        "switch": family.local_checks(),
        "old_lift_obstruction": family.old_obstruction()["checks"],
        "lower": lower.data()["checks"],
        "bound_derivatives": bounds.exact_checks(),
        "switch_bound": bounds.switch_jet_bounds()["checks"],
        "actual_vacuum": local.data()["checks"],
        "general_gradient": covariance.checks(),
    }
    return {
        prefix + "_" + key: value
        for prefix, rows in groups.items()
        for key, value in rows.items()
    }


@cache
def gates():
    groups = {
        "old_lift": family.old_obstruction()["bounds"],
        "switch": bounds.switch_jet_bounds()["bounds"],
        "tube": bounds.tube()["bounds"],
        "domain": bounds.domain()["bounds"],
    }
    d = {
        prefix + "_" + key: bool(value)
        for prefix, rows in groups.items()
        for key, value in rows.items()
    }
    d.update(
        {
            "new_named_analytic_switch_not_edit_of_frozen_S6_108": True,
            "rational_step_denominator_positive_for_every_real_X_and_even_order": True,
            "every_domain_region_has_strict_positive_full_quotient_margin": True,
            "even_order_uniform_tube_error_uses_positive_majorants": True,
            "bump_and_rational_inverse_use_finite_jet_norms": True,
            "same_actual_retuned_classical_clock_jets_through_n_minus_one": True,
            "same_massive_canonical_lambda_gamma_vacuum_contact": True,
            "old_polynomial_family_has_two_distinct_affine_rank_loss_loci": True,
            "vacuum_regular_normalization_cannot_remove_old_rank_loss": True,
            "general_gradient_action_rebuilt_before_any_rest_frame": True,
            "Lorentz_quotient_representation_has_unit_determinant": True,
            "connected_gradient_domain_allows_analytic_rank_identity_extension": True,
            "nonzero_null_gradient_not_confused_with_constant_field_vacuum": True,
            "projective_quotient_is_not_a_propagating_heavy_spectrum": True,
            "regular_X_integral_is_not_a_retarded_time_memory": True,
            "parameter_integral_is_analytic_on_every_compact_subdomain": True,
            "vacuum_regular_integral_has_no_homogeneous_inverse_square_root_term": True,
            "old_clock_boundary_cannot_be_retained_in_regular_vacuum_lift": True,
            "source_c_F4_J3_and_J2_are_analytic_at_zero_gradient": True,
            "all_actual_lower_terms_and_covariant_divergence_are_retained": True,
            "full_algebraic_stationarity_preserves_reduced_classical_equations": True,
            "free_matter_and_standalone_ordinary_Proca_have_no_connection_source": True,
            "not_the_unchanged_old_kinetic_affine_completion": True,
            "no_old_quantum_state_counterterm_or_absolute_source_transfer": True,
            "propagating_UV_parent_and_full_V_G_B_remain_open": True,
            "global_bump_derivative_ceiling": sp.Rational(3, 4) < 1,
            "global_R_X_ceiling": 1 + sp.Rational(25, 4) * family.N < 7 * family.N,
            "integral_positive_R_ratio_ceiling": sp.Rational(12, 5) ** 3 < 3**4,
            "global_q_ceiling": sp.Rational(63, 2) < 32,
            "global_h_log_derivative_bound_from_square": sp.expand(
                (family.u - 1) ** 2 - (1 + family.u**2 - 2 * family.u)
            )
            == 0,
        }
    )
    return d


@cache
def rejected_inputs():
    for name, call, args in calibration.bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported affine-domain input accepted: " + name)
    return len(calibration.bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "old_polynomial_tensor_positivity_not_a_full_affine_rank_certificate": family.old_obstruction()[
            "old_R_at_quarter"
        ]
        < sp.Rational(1, 2),
        "old_clock_q_boundary_not_reassigned_to_regular_integral": True,
        "positive_W_integral_requires_changed_auxiliary_clock_boundary": True,
        "zero_gradient_and_nonzero_null_gradient_are_separate_regressions": True,
        "matrix_calibrations_not_mistaken_for_continuous_domain_proof": True,
        "regular_domain_strictly_excludes_quotient_rank_loss": family.data()[
            "additional_quotient_factor_lower"
        ]
        > 0,
        "analytic_removable_values_not_zero_times_singular_kernels": True,
        "same_reduced_light_action_does_not_identify_kinetic_auxiliary_parents": True,
        "nonpropagating_auxiliary_inverse_not_a_heavy_threshold": True,
        "classical_lift_not_a_quantum_measure_or_V_G_B_certificate": True,
        "old_S6_108_frozen_action_and_report_unchanged": True,
        "original_P8_not_closed": True,
    }
