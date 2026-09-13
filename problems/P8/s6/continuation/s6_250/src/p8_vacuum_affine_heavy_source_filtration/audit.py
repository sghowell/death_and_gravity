"""Exact source-only boundary, unpromoted original frontier and unsupported claims."""

from functools import cache

import sympy as s
from p8_vacuum_affine_two_mass_nonscalar_inverse import audit as previous

from . import graphs, influence, source

ITEM = {
    "id": "QG2_H8A420_complete_prepared_heavy_source_influence_canonical_clock_jet_and_countergraded_source_loop_filtration",
    "status": "EXACT_CONDITIONAL_GAUSSIAN_SOURCE_INFLUENCE_AND_FORMAL_ANCESTRY_PRESERVING_LOOP_LOWER_BOUNDS_NOT_CURVED_COUNTERFUNCTIONAL_OMITTED_LOOP_ERRORS_NONLINEAR_UV_REGGE_OR_P8",
}


def require_parameters(kappa, mass_squared, cubic, localizer):
    wanted = (source.KAPPA, source.MASS2, source.G, source.A)
    for got, expected in zip(
        (kappa, mass_squared, cubic, localizer), wanted, strict=True
    ):
        if (
            isinstance(got, bool)
            or not isinstance(got, (int, s.Integer, s.Rational))
            or got != expected
        ):
            raise ValueError("Keep all four exact current physical source parameters")
    return wanted


def require_claim(label):
    if not isinstance(label, str) or label not in (
        "exact_conditional_prepared_H_source_influence",
        "generic_total_field_clock_source_filtration",
        "formal_countergraded_source_graph_lower_bounds",
        "complete_physical_forced_KG_on_unchanged_FLRW",
    ):
        raise ValueError("Keep the explicit source-only observable and formal boundary")
    return label


def require_preparation(label):
    if (
        not isinstance(label, str)
        or label != "unchanged_zero_mean_Gaussian_H_with_original_source_germ"
    ):
        raise ValueError("Keep the actual zero-mean Gaussian heavy preparation")
    return label


def require_subtraction(label):
    if (
        not isinstance(label, str)
        or label
        != "source_ancestor_and_background_jet_preserving_formal_countervertices_only"
    ):
        raise ValueError(
            "No missing curved counterfunctional or new finite term is supplied"
        )
    return label


def require_contour(label):
    if (
        not isinstance(label, str)
        or label != "complete_two_history_CTP_with_retarded_physical_difference_leg"
    ):
        raise ValueError(
            "A one-copy retarded action is not the physical influence functional"
        )
    return label


def frontier():
    return previous.frontier()


def matching():
    return [*previous.matching(), dict(ITEM)]


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("No original or previous matching obligation is promoted")
    return True


@cache
def packets():
    return {
        "whole_actual_source_and_physical_KG": source.data(),
        "entire_prepared_Gaussian_source_influence": influence.data(),
        "full_countergraded_source_graph_filtration": graphs.data(),
    }


@cache
def residuals():
    return {
        name + "_" + key: s.ImmutableMatrix(value.applyfunc(s.cancel))
        if isinstance(value, s.MatrixBase)
        else s.cancel(value)
        for name, data in packets().items()
        for key, value in data["checks"].items()
    }


def scalar_entry_count():
    return sum(
        len(value) if isinstance(value, s.MatrixBase) else 1
        for value in residuals().values()
    )


@cache
def gates():
    return {
        **{
            name + "_" + key: bool(value)
            for name, data in packets().items()
            for key, value in data["gates"].items()
        },
        "all_nine_original_primitive_rows_retained": len(frontier()) == 9,
        "all_matching_identifiers_unique": len({r["id"] for r in matching()})
        == len(matching()),
        "no_interacting_state_or_missing_counterfunctional_inferred": True,
        "source_zeros_not_a_bound_for_other_loops_or_nonlinear_feedback": True,
    }


def observable():
    return {
        "exact": "The complete current physical H source with its fixed finite extension has a canonical clock jet of degree8. Exact prepared H integration keeps the full metric determinant and both source contour kernels; the source part starts degree16 and the physical coherent source response has its full KG normalization.",
        "formal": "For generic source-dependent nonheavy graphs, E+2Ltotal>=16; a one-heavy mean first becomes topologically permitted at loop4. Countervertex ancestry preserves weighted degree under complete subgraph contraction.",
        "not_a_counterfunctional": "The structural renormalization statement is conditional on full source-ancestor/background-jet-preserving subtraction. The missing complete curved light/gravity/source counterfunctional and a full interacting initial state are not constructed.",
        "mandatory": "Source-independent heavy determinant, light/metric/M1 and mixed gravitational loops remain. Neither coefficients nor omitted-loop errors are bounded by graph counting.",
        "remaining": "No source-free quantum correction bound, evaluated coupled nonlinear feedback, same-state bounce, physical gravitational limit, UV/Regge or original P8 closure.",
    }


def bad_cases():
    invalid = (
        True,
        False,
        1.0,
        s.Float(1),
        "1",
        None,
        s.oo,
        -s.oo,
        s.zoo,
        s.I,
        s.nan,
        s.Symbol("x"),
    )
    templates = (
        (
            "parameters",
            require_parameters,
            (source.KAPPA, source.MASS2, source.G, source.A),
        ),
        ("claim", require_claim, ("exact_conditional_prepared_H_source_influence",)),
        (
            "preparation",
            require_preparation,
            ("unchanged_zero_mean_Gaussian_H_with_original_source_germ",),
        ),
        (
            "subtraction",
            require_subtraction,
            (
                "source_ancestor_and_background_jet_preserving_formal_countervertices_only",
            ),
        ),
        (
            "contour",
            require_contour,
            ("complete_two_history_CTP_with_retarded_physical_difference_leg",),
        ),
        ("external_count", graphs.minimum_loops, (2, "pure_nonheavy_source")),
    )
    out = []
    for i, value in enumerate(invalid):
        for name, call, args in templates:
            for j in range(len(args)):
                bad = list(args)
                bad[j] = value
                out.append((f"{name}_invalid_{i}_{j}", call, tuple(bad)))
    for i in range(4):
        values = [source.KAPPA, source.MASS2, source.G, source.A]
        values[i] += 1
        out.append(
            (f"changed_current_parameter_{i}", require_parameters, tuple(values))
        )
    for value in (
        "all_loops_small",
        "nonzero_renormalized_amplitude",
        "complete_curved_counterfunctional",
        "full_interacting_state",
        "tensor_first_contribution_exactly_seven",
        "source_free_gravity_loops_absent",
        "vacuum_source_order_nine",
        "nonlinear_P8_closed",
        "coupled_inverse_overwhelmed_by_exp_minus_A",
    ):
        out.append(("unsupported_" + value, require_claim, (value,)))
    for i in range(9):
        rows = frontier()
        rows[i]["status"] = "COMPLETE"
        out.append((f"primitive_promotion_{i}", validate_scope, (rows, matching())))
    for i in range(len(matching())):
        rows = matching()
        rows[i]["status"] = "COMPLETE"
        out.append((f"matching_promotion_{i}", validate_scope, (frontier(), rows)))
    out.extend(
        (
            ("missing_primitive", validate_scope, (frontier()[:-1], matching())),
            ("missing_matching", validate_scope, (frontier(), matching()[:-1])),
            ("zero_external_H_sector", graphs.minimum_loops, (0, "one_heavy_source")),
            (
                "negative_external_count",
                graphs.minimum_loops,
                (-1, "pure_nonheavy_source"),
            ),
            ("wrong_one_copy_contour", require_contour, ("half_J_retarded_J",)),
            (
                "coherent_mean_reset",
                require_preparation,
                ("arbitrary_coherent_H_state",),
            ),
        )
    )
    return out


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (ValueError, TypeError):
            continue
        raise ValueError("Unsupported source-filtration claim accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "one_copy_retarded_action_gives_advanced_contamination": influence.data()[
            "gates"
        ]["one_copy_retarded_quadratic_is_not_the_causal_force"],
        "full_source_free_heavy_determinant_retained": influence.data()["gates"][
            "source_free_full_heavy_determinant_not_removed"
        ],
        "countervertex_loop_grade_not_erased": graphs.weighted_vertex(1, 4)
        != graphs.weighted_vertex(1, 0),
        "one_heavy_mean_can_precede_pure_source_background_effect": graphs.minimum_loops(
            1, "one_heavy_source"
        )
        < graphs.minimum_loops(1, "pure_nonheavy_source"),
        "actual_source_generic_eighth_jet_nonzero": source.data()["gates"][
            "entire_order_eight_source_not_zero_function"
        ],
        "individual_spatial_projection_not_generic_nonzero_amplitude": True,
        "source_free_mixed_gravity_channel_still_required": True,
        "all_original_primitive_statuses_unchanged": frontier() == previous.frontier(),
    }
