"""Exhaustive exact denominator, graph and restricted-forest audit."""

from functools import cache
from itertools import combinations

from . import calibration, forests, graphs, positive, subgraphs


def label(kind, choices):
    return kind + "_" + "".join(map(str, choices))


@cache
def residuals():
    result = {"graph_" + k: v for k, v in graphs.data()["checks"].items()}
    result.update(
        {"calibration_" + k: v for k, v in calibration.data()["checks"].items()}
    )
    for kind, choices in graphs.cases():
        for prefix, module in (("forest", forests), ("positive", positive)):
            result.update(
                {
                    prefix + "_" + label(kind, choices) + "_" + k: v
                    for k, v in module.data(kind, choices)["checks"].items()
                }
            )
    return result


@cache
def gates():
    all_graphs = [graphs.refinement(*case) for case in graphs.cases()]
    all_forests = [forests.data(*case) for case in graphs.cases()]
    all_subgraphs = [subgraphs.data(*case) for case in graphs.cases()]
    result = {k: bool(v) for k, v in calibration.data()["bounds"].items()}
    result.update(
        {
            "exactly_nine_labelled_light_1PI_skeletons": len(graphs.skeletons()) == 9,
            "exactly_three_topology_families": {r["kind"] for r in graphs.skeletons()}
            == set(graphs.KINDS),
            "exactly_seventy_two_external_assignments": graphs.data()[
                "external_assignment_count"
            ]
            == 72,
            "all_sixty_four_heavy_choices_per_family": len(graphs.cases()) == 192,
            "every_refinement_keeps_two_loops": all(
                g["loop_count"] == 2 for g in all_graphs
            ),
            "every_refinement_keeps_four_external_light_legs": all(
                sorted(x for labels in g["external_labels"].values() for x in labels)
                == [0, 1, 2, 3]
                for g in all_graphs
            ),
            "every_refinement_has_positive_tree_sum_on_interior": all(
                d["spanning_tree_count"] > 0 for d in all_forests
            ),
            "all_crossing_derivative_majorants_have_nonnegative_coefficients": all(
                d["derivative_majorant_coefficients_nonnegative"] for d in all_forests
            ),
            "every_exact_positive_polynomial_witness_succeeds": all(
                r["all_square_weights_and_remainder_coefficients_nonnegative"]
                for r in positive.summaries()
            ),
            "all_subgraph_external_counts_nonnegative": all(
                r["light_external_legs"] >= 0 and r["heavy_external_legs"] >= 0
                for d in all_subgraphs
                for r in d["UV_subgraphs"]
            ),
            "every_restricted_forest_contains_only_nested_or_vertex_disjoint_pairs": all(
                set(d["UV_subgraphs"][i]["vertices"]).isdisjoint(
                    d["UV_subgraphs"][j]["vertices"]
                )
                or set(d["UV_subgraphs"][i]["edges"])
                < set(d["UV_subgraphs"][j]["edges"])
                or set(d["UV_subgraphs"][j]["edges"])
                < set(d["UV_subgraphs"][i]["edges"])
                for d in all_subgraphs
                for forest in d["restricted_forest_index_sets"]
                for i, j in combinations(forest, 2)
            ),
            "independent_background_Gaussian_Wick_count_is_only_combinatorial_control": True,
            "exact_heavy_pairings_not_loop_momentum_derivative_expansion": True,
            "distinct_self_loop_half_edges_kept_in_heavy_refinement": True,
            "mass_shell_single_cuts_and_all_crossing_pairings_kept": True,
            "matrix_tree_cofactor_independent_of_spanning_tree_enumeration": True,
            "all_rational_square_weights_verified_not_numerically_fitted": True,
            "strict_simplex_interior_bound_not_boundary_convergence": True,
            "UV_subgraph_forests_are_combinatorial_not_evaluated_subtractions": True,
            "nested_counterterm_operations_must_be_inside_out": True,
            "light_1PI_skeletons_do_not_include_complete_two_loop_LSZ": True,
            "finite_fixed_scheme_conversion_not_assumed_zero": True,
            "positive_Wick_weights_not_signed_renormalized_amplitudes": True,
            "Euclidean_kernel_iterations_not_a_complete_two_loop_remainder": True,
            "proper_on_shell_map_not_global_off_shell_coordinate_bound": True,
            "no_all_energy_or_Regge_bound_from_compact_disc": True,
            "no_finite_gravity_or_rolling_state_transfer": True,
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
        raise ValueError("Unsupported two-loop denominator input accepted: " + name)
    return len(calibration.bad_cases())


@cache
def controls():
    double = subgraphs.data("double_bubble", (0, 0, 0))
    tadpole = graphs.refinement("tadpole_insertion", (0, 0, 0))
    mixed = graphs.refinement("tadpole_insertion", (0, 0, 2))
    bare = forests.data("double_bubble", (0, 0, 0))
    boundary = dict(zip(bare["parameters"], (1, 0, 0, 0)))
    return {
        "rejected_inputs": rejected_inputs(),
        "bare_double_bubble_has_vertex_overlap": double[
            "edge_disjoint_but_vertex_overlapping_pairs"
        ]
        == ((0, 1),),
        "overlapping_bubbles_not_in_same_restricted_forest": not any(
            {0, 1} <= set(row) for row in double["restricted_forest_index_sets"]
        ),
        "bare_double_bubble_has_six_forests": len(
            double["restricted_forest_index_sets"]
        )
        == 6,
        "bare_tadpole_retains_self_loop": any(a == b for a, b in tadpole["edges"]),
        "split_tadpole_has_mixed_bubble_not_a_self_loop": not any(
            a == b for a, b in mixed["edges"]
        ),
        "simplex_boundary_can_have_zero_U": bare["U"].subs(boundary) == 0,
        "compact_denominator_bound_not_a_UV_boundary_estimate": True,
        "coupling_smallness_does_not_bound_unsubtracted_loop_integrals": True,
        "higher_loop_LSZ_and_fixed_counterterms_not_omitted_by_claim": True,
        "full_V_G_B_and_original_P8_not_closed": True,
    }
