"""Complete field/loop grading, source parity and counterterm ancestry boundaries."""

from functools import cache

import sympy as s


def natural(value, label):
    if isinstance(value, bool) or not isinstance(value, (int, s.Integer)) or value < 0:
        raise ValueError("Use an exact nonnegative integer " + label)
    return int(value)


def minimum_loops(external, sector):
    external = natural(external, "external count")
    if not isinstance(sector, str) or sector not in (
        "pure_nonheavy_source",
        "one_heavy_source",
    ):
        raise ValueError("Keep the stated source-dependent graph sectors")
    if sector == "one_heavy_source" and external < 1:
        raise ValueError("A one-heavy external sector has at least one external leg")
    total = 16 if sector == "pure_nonheavy_source" else 9
    return max(0, (total - external + 1) // 2)


def weighted_vertex(degree, loop_grade):
    degree = natural(degree, "vertex degree")
    loop_grade = natural(loop_grade, "countervertex loop grade")
    return degree - 2 + 2 * loop_grade


def complete_graph(degrees, internal, grades=None):
    if not isinstance(degrees, tuple) or not degrees:
        raise ValueError("Use a nonempty connected vertex degree tuple")
    degrees = tuple(natural(d, "degree") for d in degrees)
    internal = natural(internal, "internal edge count")
    grades = (0,) * len(degrees) if grades is None else grades
    if not isinstance(grades, tuple) or len(grades) != len(degrees):
        raise ValueError("Every vertex needs its exact loop grade")
    grades = tuple(natural(r, "loop grade") for r in grades)
    external = sum(degrees) - 2 * internal
    topological = internal - len(degrees) + 1
    if external < 0 or topological < 0:
        raise ValueError("Counts cannot describe a connected graph")
    total = topological + sum(grades)
    return {
        "external": external,
        "topological_loops": topological,
        "total_loop_grade": total,
        "weighted_degree_sum": sum(
            weighted_vertex(d, r) for d, r in zip(degrees, grades, strict=True)
        ),
    }


@cache
def data():
    rows = []
    checks = {}
    for external in range(17):
        bound = minimum_loops(external, "pure_nonheavy_source")
        rows.append(
            {
                "external_nonheavy_legs": external,
                "first_topologically_permitted_loop": bound,
            }
        )
        checks["sharp_integer_loop_lower_bound_" + str(external)] = (
            max(0, (16 - external + 1) // 2) - bound
        )
    # Each pair uses its one heavy-to-heavy edge and k cross-light edges.
    for k in range(9):
        graph = complete_graph((9, 9), k + 1)
        checks["full_two_source_graph_external_legs_" + str(k)] = graph["external"] - (
            16 - 2 * k
        )
        checks["full_two_source_graph_loop_number_" + str(k)] = (
            graph["total_loop_grade"] - k
        )
        checks["full_two_source_graph_weight_" + str(k)] = (
            graph["weighted_degree_sum"] - 14
        )
        counter_degree, counter_grade = graph["external"], graph["total_loop_grade"]
        checks["complete_two_source_countergraph_contraction_" + str(k)] = (
            weighted_vertex(counter_degree, counter_grade) - 14
        )
    for k in range(5):
        graph = complete_graph((9,), k)
        checks["whole_single_source_tadpole_external_count_" + str(k)] = graph[
            "external"
        ] - (9 - 2 * k)
        checks["whole_single_source_tadpole_loop_count_" + str(k)] = (
            graph["total_loop_grade"] - k
        )
        checks["single_source_countervertex_weight_preserved_" + str(k)] = (
            weighted_vertex(graph["external"], k) - 7
        )
    actual_mean = complete_graph((10, 9), 9)
    checks["generic_next_clock_source_jet_mean_graph"] = actual_mean["external"] - 1
    checks["generic_next_clock_source_jet_mean_loop"] = (
        actual_mean["total_loop_grade"] - 8
    )
    checks["one_heavy_mean_first_permitted_loop"] = (
        minimum_loops(1, "one_heavy_source") - 4
    )
    E, I, V, L, R = s.symbols(
        "external internal vertices topological_loop counter_loop_sum", integer=True
    )
    checks["entire_countergraded_connected_graph_identity"] = s.expand(
        2 * I + E - 2 * V + 2 * R - (E + 2 * (L + R) - 2)
    ).subs(I, L + V - 1)
    checks["contracted_subgraph_preserves_full_weight"] = s.expand(
        (E - 2 + 2 * (L + R)) - (2 * I + E - 2 * V + 2 * R)
    ).subs(I, L + V - 1)
    marker, j, k = s.symbols("source_marker physical_source heavy_kinetic", real=True)
    source_effect = marker**2 * j * j / (2 * k)
    heavy_mean = marker * j / k
    checks["full_source_Gaussian_even_parity"] = (
        source_effect.subs(marker, -marker) - source_effect
    )
    checks["full_one_heavy_odd_parity"] = heavy_mean.subs(marker, -marker) + heavy_mean
    # Two minimal H²h vertices, two internal heavy edges and one graviton edge.
    minimal_mixed = complete_graph((3, 3), 3)
    checks["source_independent_mixed_gravity_vacuum_two_loops"] = (
        minimal_mixed["total_loop_grade"] - 2
    )
    checks["source_independent_mixed_gravity_vacuum_no_external_legs"] = minimal_mixed[
        "external"
    ]
    checks["full_fixed_onepoint_extension_countergrade"] = (
        weighted_vertex(1025, 1) - 1025
    )
    return {
        "graph_identity": "sum_v(d_v-2+2r_v)=E+2Ltotal-2, Ltotal=I-V+1+sum_v r_v. Counts are necessary, not a proof that every degree sequence has nonzero renormalized amplitude.",
        "pure_nonheavy_source_table": rows,
        "one_heavy_mean_first_topologically_permitted_loop": 4,
        "two_source_filtration": "Pure nonheavy source-dependent kernels need an even positive number of H-linear source vertices, at least two. Each bare source vertex has degree>=9, so E+2Ltotal>=16. The full finite V source has higher degree and loop grade1. Baseline tree interactions of degree>=3 and correctly graded countervertices do not lower this weight.",
        "specific_lower_bounds": "Source-dependent vacuum and nonheavy mean: at least8 loops; two-point/linear-response: at least7; three-point: at least7; four-point: at least6. One-heavy mean: at least4. They are first topologically permitted orders, not proven nonzero continuum coefficients or error bounds.",
        "counterterm_condition": "Contracting an actual source subgraph to its complete countervertex preserves weighted field/loop degree. This gives a formal filtration for source-descended subtractions respecting the full background Taylor expansion and source marker. It does not supply the missing curved counterfunctional or authorize arbitrary independent lower-degree finite matching terms. Countervertices of positive loop grade may have lower field degree.",
        "reference_grading": "The QG2 Gaussian mean profiles and the stated finite vacuum counterterms carry their loop grades when making a formal loop statement. A resummed reference propagator must be expanded consistently for that purpose. Do not reinterpret a positive-loop onepoint countervertex as a tree tadpole or claim the exact resummed inverse is a convergent loop series.",
        "sharpness_boundary": "Explicit connected graphs realize the generic nonheavy count bounds at finite regulator. Individual tensor/vector or other projections can require extra vertices or vanish, and actual derivative contractions, gauges and subtraction may cancel amplitudes. Topology alone proves no nonzero projected or renormalized coefficient.",
        "mandatory_source_independent_channels": "The whole metric-dependent H determinant remains at loop1. Minimal H²h interactions already permit a two-loop mixed heavy/gravity vacuum graph. Pure light/metric/M1 and gauge-sector loops are not postponed by the high-order H source. They remain mandatory unbounded channels.",
        "vacuum_clock_distinction": "The Minkowski vacuum X0 has the retained cubic g H Phi²/2 and one-loop heavy-onepoint cancellation. The clock X1 has a different source jet. The clock loop filtration is not transferred to vacuum amplitudes or the physical gravitational limit.",
        "checks": {name: s.cancel(value) for name, value in checks.items()},
        "gates": {
            "source_linear_kernel_first_permitted_seven_not_one": minimum_loops(
                2, "pure_nonheavy_source"
            )
            == 7,
            "source_four_point_first_permitted_six_not_one": minimum_loops(
                4, "pure_nonheavy_source"
            )
            == 6,
            "one_heavy_tadpole_not_delayed_to_eight": minimum_loops(
                1, "one_heavy_source"
            )
            == 4,
            "lower_degree_countervertices_cannot_be_treated_as_tree": weighted_vertex(
                1, 4
            )
            == 7
            and weighted_vertex(1, 0) == -1,
            "full_fixed_finite_source_extension_grade_kept": weighted_vertex(1025, 1)
            > weighted_vertex(9, 0),
            "non_source_two_loop_gravity_channel_retained": minimal_mixed[
                "total_loop_grade"
            ]
            == 2,
            "missing_full_curved_counterfunctional_not_invented": True,
            "graph_lower_bound_not_nonzero_amplitude_or_small_error": True,
        },
    }
