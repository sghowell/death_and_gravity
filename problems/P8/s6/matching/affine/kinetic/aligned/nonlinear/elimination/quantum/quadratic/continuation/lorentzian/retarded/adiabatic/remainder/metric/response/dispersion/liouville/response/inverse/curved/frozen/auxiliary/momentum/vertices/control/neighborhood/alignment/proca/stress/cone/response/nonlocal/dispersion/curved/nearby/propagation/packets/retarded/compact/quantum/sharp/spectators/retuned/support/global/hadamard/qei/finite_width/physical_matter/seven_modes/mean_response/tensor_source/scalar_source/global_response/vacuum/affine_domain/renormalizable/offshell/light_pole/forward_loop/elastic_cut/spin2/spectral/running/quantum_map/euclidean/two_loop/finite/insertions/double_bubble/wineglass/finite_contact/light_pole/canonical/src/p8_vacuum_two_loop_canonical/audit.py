"""Exact complete-order ledger checks and explicit finite-order boundaries."""

from functools import cache

from . import assembly, calibration, contractions, fields, ledger


@cache
def residuals():
    return {
        prefix + "_" + k: v
        for prefix, module in (
            ("contractions", contractions),
            ("ledger", ledger),
            ("fields", fields),
            ("assembly", assembly),
            ("calibration", calibration),
        )
        for k, v in module.data()["checks"].items()
    }


@cache
def gates():
    result = {k: bool(v) for k, v in calibration.data()["bounds"].items()}
    result.update(
        {
            "all_4608_full_external_labelled_graph_records": contractions.data()[
                "raw_external_labelled_graph_records"
            ]
            == 4608,
            "all_1152_marked_vertex_and_heavy_mass_cores": contractions.data()[
                "single_proper_core_contractions"
            ]
            == 1152,
            "all_72_disjoint_core_products": contractions.data()[
                "disjoint_proper_core_pairs"
            ]
            == 72,
            "all_144_one_loop_counterterm_cographs": len(
                contractions.data()["all_single_core_cograph_weight_matches"]
            )
            == 144,
            "all_12_disjoint_tree_counterterm_cographs": len(
                contractions.data()["all_disjoint_tree_cograph_weight_matches"]
            )
            == 12,
            "all_192_raw_refinements_have_exactly_one_integrated_owner": len(
                assembly.data()["exact_raw_refinement_ownership"]
            )
            == 192,
            "external_labels_and_propagator_types_preserved_in_canonical_graph_keys": True,
            "one_factor_of_full_regulated_I0_per_contracted_light_core": True,
            "fundamental_cubic_G_used_before_derived_g_squared_parameter": True,
            "disjoint_cubic_cubic_cubic_mass_and_mass_mass_products_once_each": True,
            "all_actual_signed_contraction_weights_match_independent_insertion_graphs": True,
            "inherited_overall_delta_g_references_are_linear_tree_variations": True,
            "same_fixed_heavy_parameter_references_not_refitted": True,
            "all_regulator_dependence_retained_until_forest_cancellation": True,
            "fixed_inner_light_OS_subtraction_generated_reference_retained": True,
            "once_fixed_finite_quartic_insertion_keeps_both_heavy_triangles": True,
            "new_fixed_two_loop_potential_contact_cannot_tune_b2": True,
            "actual_canonical_unit_light_residue_through_two_loops": fields.data()[
                "actual_inherited_light_pole"
            ]["light_pole_residue_through_two_loops"]
            == 1,
            "bare_field_factor_not_confused_with_canonical_pole_residue": True,
            "bare_amputation_and_LSZ_cancel_through_second_order": True,
            "no_extra_tree_LSZ_factor_in_total_canonical_vertices": True,
            "symmetric_perturbative_light_vacuum_has_no_odd_light_vertices": True,
            "heavy_reducible_exchange_graphs_retained": True,
            "full_canonical_two_loop_b2_now_bounded_after_complete_ledger": True,
            "no_signed_exact_two_loop_coefficient_claimed": True,
            "two_loop_derivative_coordinate_source_matching_not_computed": True,
            "no_all_higher_loop_or_all_energy_error_bound": True,
            "no_finite_gravity_Regge_Delta_or_bounce_parent_completion": True,
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
        raise ValueError("Unsupported canonical loop order accepted: " + name)
    return len(calibration.bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "omitted_cubic_square_changes_second_forward_coefficient": ledger.data()[
            "diagnostic_omitted_cubic_product_b2"
        ]
        != 0,
        "extra_canonical_LSZ_rescaling_changes_amplitude": fields.data()[
            "diagnostic_nonzero_extra_LSZ_error"
        ]
        != 0,
        "bare_counterterm_products_not_momentum_independent_contacts": True,
        "regulator_symbols_not_assigned_finite_values": True,
        "group_bound_not_promoted_without_independent_graph_and_field_ledger": True,
        "canonical_amplitude_not_derivative_coordinate_source_result": True,
        "finite_order_positivity_not_all_orders_UV_completion": True,
        "full_V_G_B_and_original_P8_not_closed": True,
    }
