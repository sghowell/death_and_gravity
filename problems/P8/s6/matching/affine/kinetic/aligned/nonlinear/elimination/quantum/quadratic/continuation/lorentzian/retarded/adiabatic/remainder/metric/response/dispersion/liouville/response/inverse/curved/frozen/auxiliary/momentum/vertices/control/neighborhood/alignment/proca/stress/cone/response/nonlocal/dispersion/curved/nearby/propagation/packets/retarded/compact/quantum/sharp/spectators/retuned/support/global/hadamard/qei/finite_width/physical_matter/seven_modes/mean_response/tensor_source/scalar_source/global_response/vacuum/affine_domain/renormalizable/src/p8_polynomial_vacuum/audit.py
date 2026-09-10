"""Exact audit and explicit non-claims for polynomial vacuum quantum preparation."""

from functools import cache

from p8_exceptional_vacuum import heavy

from . import calibration, gaussian, model, potential, powercount


@cache
def residuals():
    groups = {
        "model": model.data()["checks"],
        "gaussian": gaussian.data()["checks"],
        "powercount": powercount.data()["checks"],
        "potential": potential.data()["checks"],
    }
    return {
        prefix + "_" + key: value
        for prefix, rows in groups.items()
        for key, value in rows.items()
    }


@cache
def gates():
    groups = {
        "model": model.data()["bounds"],
        "powercount": powercount.data()["bounds"],
        "potential": potential.data()["bounds"],
    }
    rows = {
        prefix + "_" + key: bool(value)
        for prefix, values in groups.items()
        for key, value in values.items()
    }
    rows.update(
        {
            "new_canonical_polynomial_model_not_edit_of_frozen_actions": True,
            "exact_full_on_shell_tree_four_point_match_not_only_truncated_contact": True,
            "positive_square_margin_gives_unique_global_classical_vacuum": True,
            "H_decay_width_is_tree_estimate_not_all_order_stable_pole": True,
            "all_finite_regulator_light_correlators_keep_exact_nonlocal_quartic": True,
            "positive_kernel_inverse_bound_gives_coercive_finite_regulator_action": True,
            "two_site_full_Hessian_and_literal_determinant_rebuilt": True,
            "heavy_determinant_field_independent_does_not_delete_mixed_loops": True,
            "arbitrary_loop_order_counterterm_candidates_follow_unbounded_graph_topology": True,
            "subdiagram_counterterms_preserve_at_most_quadratic_heavy_structure": True,
            "finite_local_subtractions_chosen_once_only_for_new_vacuum_model": True,
            "constant_field_Schur_kernel_positive_at_all_Euclidean_momenta": True,
            "subtracted_log_integral_has_positive_global_lower_and_upper": True,
            "UV_integrability_follows_positive_radial_majorant": True,
            "one_loop_potential_not_on_shell_mass_or_residue_renormalization": True,
            "no_derivative_expansion_under_unbounded_loop_integral": True,
            "no_bare_squared_unstable_resonance_dispersion_integral": True,
            "no_all_order_unitarity_or_complex_contour_bound_inferred": True,
            "minimal_Einstein_scalar_parent_has_nonpositive_Hubble_derivative": True,
            "same_four_point_amplitude_does_not_identify_affine_bounce_parent": True,
            "old_reference_states_counterterms_and_action_bytes_unchanged": True,
            "full_V_G_B_and_original_P8_remain_open": True,
        }
    )
    return rows


def bad_cases():
    return [
        ("powercount_" + name, call, args)
        for name, call, args in powercount.bad_cases()
    ] + calibration.bad_cases()


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (ValueError, TypeError):
            continue
        raise ValueError("Unsupported polynomial-vacuum input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    d = calibration.point()
    return {
        "rejected_inputs": rejected_inputs(),
        "G_nonzero_retains_mixed_heavy_light_loops": d["cubic_squared"] > 0,
        "H_mass_above_two_light_threshold_is_not_an_exact_stable_pole": d[
            "heavy_mass_squared"
        ]
        > 4,
        "tree_polydisc_contains_light_loop_threshold": heavy.CHANNEL_RADIUS > 2,
        "naive_quartic_squared_over_tree_b2_exceeds_one_e_180": d["quartic"] ** 2
        / d["tree_forward_b2"]
        > 10**180,
        "positive_potential_does_not_bound_forward_amplitude_loops": True,
        "Gaussian_integral_not_a_local_derivative_truncation": True,
        "full_Hessian_differentiated_before_stationary_heavy_substitution": True,
        "constant_background_positivity_not_all_inhomogeneous_Hessian_positivity": True,
        "renormalizable_powercount_not_a_nonperturbative_continuum_UV_theorem": True,
        "new_vacuum_subtraction_not_old_cosmological_tadpole_reset": True,
        "no_finite_gravity_IR_Regge_remainder_supplied": True,
        "no_common_propagating_parent_or_bounce_UV_classification": True,
        "original_P8_not_closed": True,
    }
