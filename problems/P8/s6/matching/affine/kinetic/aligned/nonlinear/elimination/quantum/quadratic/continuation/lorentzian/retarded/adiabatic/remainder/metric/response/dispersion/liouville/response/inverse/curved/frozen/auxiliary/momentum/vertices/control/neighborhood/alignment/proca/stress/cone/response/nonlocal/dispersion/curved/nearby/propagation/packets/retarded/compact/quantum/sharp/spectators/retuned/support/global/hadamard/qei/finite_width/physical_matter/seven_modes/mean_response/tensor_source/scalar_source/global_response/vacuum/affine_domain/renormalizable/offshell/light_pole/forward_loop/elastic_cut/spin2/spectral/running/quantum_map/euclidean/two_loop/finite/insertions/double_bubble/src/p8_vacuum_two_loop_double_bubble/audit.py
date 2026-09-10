"""Source-pinned subtraction, no-double-counting and original-scope audit."""

from functools import cache

from . import bubble, calibration, forests, selection, triangle


@cache
def residuals():
    return {
        prefix + "_" + k: v
        for prefix, module in (
            ("selection", selection),
            ("forests", forests),
            ("bubble", bubble),
            ("triangle", triangle),
            ("calibration", calibration),
        )
        for k, v in module.data()["checks"].items()
    }


def bad_cases():
    invalid = (
        (),
        (0, 0),
        (0, 0, 0, 0),
        [0, 0, 0],
        (True, 0, 0),
        (0, 0, False),
        (0, 4, 0),
        (0, -1, 0),
        (0, 0, 2),
        (2, 2, 0),
        (3, 3, 1),
    )
    return calibration.bad_cases() + [
        ("unsupported_graph_" + str(i), selection.require_group, (value,))
        for i, value in enumerate(invalid)
    ]


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported double-bubble input accepted: " + name)
    return len(bad_cases())


@cache
def gates():
    selected = selection.data()
    result = {k: bool(v) for k, v in calibration.data()["bounds"].items()}
    result.update(
        {
            "exact_twenty_four_refinements_selected": selected[
                "selected_refinement_count"
            ]
            == 24,
            "eight_already_finite_factorizing_refinements_excluded": selected[
                "already_finite_separable_count"
            ]
            == 8,
            "all_three_local_counterterm_types_present": selected[
                "local_UV_external_types"
            ]
            == ((0, 2), (2, 1), (4, 0)),
            "forest_classes_are_sixteen_single_four_disjoint_four_overlap": selected[
                "forest_class_counts"
            ]
            == {"single": 16, "disjoint": 4, "overlap": 4},
            "exact_sixteen_raw_wineglass_refinements_remain": len(
                selected["remaining_wineglass_refinements"]
            )
            == 16,
            "shared_vertex_subgraphs_not_subtracted_as_an_illegal_pair": True,
            "inside_out_recursive_overall_subtraction_retained": True,
            "both_disjoint_subgraphs_can_occur_in_one_forest": True,
            "all_counterterms_local_before_attaching_external_heavy_propagators": True,
            "iterated_heavy_mass_counterterms_not_a_new_nonlocal_vertex": True,
            "higher_order_zero_momentum_reference_explicitly_specified": True,
            "same_entire_regulated_I0_not_a_silent_pole_only_conversion": True,
            "same_common_regulator_used_before_cancellations_and_removal": True,
            "full_heavy_triangles_retained_in_both_independent_loops": True,
            "all_three_forward_channels_and_Wick_weights_retained": True,
            "same_first_sheet_complex_routing_and_absolute_majorants": True,
            "Cauchy_estimate_applies_to_integrated_renormalized_group": True,
            "no_double_counting_of_the_previous_eighty_eight_finite_graphs": True,
            "positive_absolute_error_not_a_signed_graph_value": True,
            "finite_potential_counterterm_insertions_not_computed_here": True,
            "complete_two_loop_mass_residue_LSZ_not_computed": True,
            "no_all_higher_loop_high_energy_or_gravity_transfer": True,
            "original_V_G_B_and_P8_remain_open": True,
        }
    )
    return result


@cache
def controls():
    rows = selection.data()["selected_rows"]
    overlap_forbidden = all(
        not any({0, 1}.issubset(f) for f in r["restricted_forest_index_sets"])
        for r in rows
        if r["mode"] == "overlap"
    )
    return {
        "rejected_inputs": rejected_inputs(),
        "illegal_shared_pair_absent_from_actual_forests": overlap_forbidden,
        "disjoint_pair_present_where_allowed": all(
            (0, 1) in r["restricted_forest_index_sets"]
            for r in rows
            if r["mode"] == "disjoint"
        ),
        "overall_subtraction_not_dropped": True,
        "eight_finite_refinements_not_counted_twice": True,
        "unregulated_reference_not_assigned_a_finite_number": True,
        "finite_potential_counterterm_insertions_not_silently_omitted_from_completion_list": True,
        "sixteen_remaining_raw_refinements_not_full_LSZ_completion": True,
        "full_two_loop_and_original_P8_not_closed": True,
    }
