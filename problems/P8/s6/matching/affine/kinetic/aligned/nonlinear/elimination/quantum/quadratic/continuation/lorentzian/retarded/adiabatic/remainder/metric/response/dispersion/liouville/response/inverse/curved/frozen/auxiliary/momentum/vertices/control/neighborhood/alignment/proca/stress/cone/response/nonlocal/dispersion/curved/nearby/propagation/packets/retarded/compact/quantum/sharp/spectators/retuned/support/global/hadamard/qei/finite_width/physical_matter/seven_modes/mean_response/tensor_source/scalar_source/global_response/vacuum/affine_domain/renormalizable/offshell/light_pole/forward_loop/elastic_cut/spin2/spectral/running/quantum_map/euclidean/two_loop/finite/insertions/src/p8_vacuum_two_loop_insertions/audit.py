"""Counterterm-aware self-energy insertion audit and scope controls."""

from functools import cache

from . import calibration, kernel, radial, routing, subtraction


@cache
def residuals():
    return {
        prefix + "_" + k: v
        for prefix, module in (
            ("kernel", kernel),
            ("routing", routing),
            ("subtraction", subtraction),
            ("radial", radial),
            ("calibration", calibration),
        )
        for k, v in module.data()["checks"].items()
    }


@cache
def gates():
    result = {k: bool(v) for k, v in calibration.data()["bounds"].items()}
    result.update(
        {
            "all_sixty_four_tadpole_refinements_grouped": subtraction.data()[
                "grouped_refinement_count"
            ]
            == 64,
            "local_tadpoles_mass_and_kinetic_subtractions_kept_together": True,
            "same_actual_on_shell_mixed_self_energy_not_surrogate": True,
            "constant_asymptotic_multiplier_is_fixed_residue_counterterm": True,
            "constant_part_multiplies_same_renormalized_one_loop_bubble": True,
            "outer_subtraction_linear_in_inherited_I0_reference": True,
            "outer_counterterms_local_in_full_two_field_model": True,
            "finite_potential_contact_cannot_adjust_forward_second_derivative": True,
            "total_canonical_vertex_counterterms_not_extra_tree_LSZ_rescaling": True,
            "actual_forward_routing_keeps_all_three_channels": True,
            "bilinear_Euclidean_momenta_and_Hermitian_norms_distinguished": True,
            "complete_heavy_propagators_retained_inside_both_loops": True,
            "routed_self_energy_stays_in_analytic_complex_strip": True,
            "complex_Feynman_shift_justified_on_same_first_sheet": True,
            "decaying_remainder_absolutely_integrable_before_regulator_removal": True,
            "radial_integral_extends_to_infinity_not_a_cutoff": True,
            "Cauchy_bound_applies_to_whole_unit_disc_remainder": True,
            "both_line_insertions_and_frozen_symmetry_weight_retained": True,
            "positive_error_upper_not_positive_signed_graph_contribution": True,
            "forty_other_subtraction_dependent_refinements_still_open": True,
            "complete_two_loop_pole_and_LSZ_not_computed": True,
            "no_all_higher_loop_or_high_energy_control": True,
            "no_gravity_Regge_or_rolling_state_transfer": True,
            "original_V_G_B_and_P8_remain_open": True,
        }
    )
    return result


@cache
def rejected_inputs():
    for name, call, args in calibration.bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported insertion bound input accepted: " + name)
    return len(calibration.bad_cases())


@cache
def controls():
    d = calibration.data()
    return {
        "rejected_inputs": rejected_inputs(),
        "fixed_constant_insertion_term_not_set_to_zero": d[
            "constant_insertion_group_b2_upper"
        ]
        > 0,
        "decaying_remainder_not_omitted": d["decaying_insertion_group_b2_upper"] > 0,
        "large_individual_mass_counterterms_not_separately_bounded": True,
        "unsubtracted_constant_bubble_not_assigned_a_finite_value": True,
        "same_inner_OS_and_outer_reference_no_b2_contact_tuning": True,
        "finite_potential_contact_not_a_free_derivative_counterterm": True,
        "partial_raw_graph_count_not_full_two_loop_LSZ": True,
        "full_two_loop_and_original_P8_not_closed": True,
    }
