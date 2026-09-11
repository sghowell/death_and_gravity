"""Connected loop grading, forest contraction and dimensional reference data."""

from functools import cache

import sympy as s

from . import profile


def new_vertex_floor(external, index=1):
    if type(external) is not int or not 0 <= external <= 10:
        raise TypeError("Require a native source-endpoint count from zero through ten")
    if type(index) is not int or index < 1:
        raise TypeError("Require a positive native new-vertex index")
    return int(s.ceiling(1 - s.Rational(external, 2) + s.Rational(8 * index + 1, 2)))


def grade(external, valences, counterterm_orders):
    if type(external) is not int or external < 0:
        raise TypeError("Require a nonnegative native endpoint count")
    if (
        not isinstance(valences, tuple)
        or not isinstance(counterterm_orders, tuple)
        or len(valences) != len(counterterm_orders)
        or not valences
    ):
        raise TypeError(
            "Require equal nonempty native valence and counterterm-order tuples"
        )
    if any(type(v) is not int or v < 0 for v in valences + counterterm_orders):
        raise TypeError("Require nonnegative native valences and loop orders")
    if (sum(valences) - external) % 2 or sum(valences) < external:
        raise ValueError("The halfedge count cannot describe a diagram")
    value = (
        1
        - s.Rational(external, 2)
        + sum(s.Rational(v - 2, 2) + j for v, j in zip(valences, counterterm_orders))
    )
    if value < 0 or value.is_integer is not True:
        raise ValueError("The effective loop grade is not a nonnegative integer")
    return int(value)


@cache
def data():
    k, eps = s.symbols("k epsilon", real=True)
    I, V, E, J = s.symbols("I V E J", real=True)
    original_weight = (2 * I + E - 2 * V) / 2 + J
    contracted_weight = (E - 2) / 2 + (I - V + 1 + J)
    t, z, r = s.symbols("t z cubic_R", real=True)
    f9 = (t * z + t**3 * r) ** 9
    dimensions = (4 - 2 * eps) - (
        (8 * k + 1) * (1 - eps) + 2 * (s.Rational(3, 2) - eps)
    )
    floors = {E: new_vertex_floor(E) for E in (0, 1, 2, 4, 8, 10)}
    return {
        "first_new_vertex": {
            "Phi_legs": 9,
            "fermion_legs": 2,
            "total_valence": 11,
            "tree_counterterm_order": 0,
            "weighted_excess": s.Rational(9, 2),
        },
        "integer_loop_lower_bounds_by_total_source_count": floors,
        "effective_loop_identity": "L_eff=1-E/2+sum_v[(n_v-2)/2+j_v]+p, where p>=0 is the extra cubic-map source/action degree when using the physical pullback. All ordinary connected nonvacuum vertices have nonnegative weighted excess; the first linear H reference counterterm has excess -1/2+1=1/2.",
        "regulator_prescription": "In D=4-2epsilon use y_D=mu^epsilon y and R_D=mu^-epsilon R. The k-th new coefficient y_D/R_D^(8k)=mu^((8k+1)epsilon)y/R^(8k). Keep this same lift before all pole subtraction and finite products.",
        "lower_bound_not_sharp_first_b2": "The table bounds the first possible changed graph, not the first nonzero b2. For a single H source the bound five is sufficient here and is not asserted sharp. Amplitudes with ten scalar fields can change already at one loop.",
        "checks": {
            "all_new_vertices_have_nine_or_more_scalar_legs": (8 * k + 1) - 1 - 8 * k,
            "first_new_total_valence": 9 + 2 - 11,
            "first_new_weighted_excess": s.Rational(11 - 2, 2) - s.Rational(9, 2),
            "four_source_floor": floors[4] - 4,
            "two_source_floor": floors[2] - 5,
            "vacuum_floor": floors[0] - 6,
            "single_source_floor_not_assumed_sharp": floors[1] - 5,
            "eight_source_lower_bound": floors[8] - 2,
            "ten_source_one_loop_change_allowed": floors[10] - 1,
            "forest_contraction_preserves_weight": s.expand(
                original_weight - contracted_weight
            ),
            "dimensional_new_coefficient": s.expand(
                dimensions - (-8 * k + (8 * k + 1) * eps)
            ),
            "first_profile_change_after_cubic_map": s.expand(f9).coeff(t, 9) - z**9,
            "no_lower_source_map_field_degrees": sum(
                s.expand(f9).coeff(t, j) for j in range(9)
            ),
            "linear_H_counterterm_nonnegative_weight": s.Rational(1 - 2, 2)
            + 1
            - s.Rational(1, 2),
            "binomial_profile_first_new_coefficient": profile.coefficient(1)
            + s.Rational(1, 8),
        },
    }
