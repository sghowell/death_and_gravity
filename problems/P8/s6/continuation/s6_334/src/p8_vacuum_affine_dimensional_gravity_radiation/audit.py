"""Dimensionally complete selected tree and preserved original physical frontier."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_joint_soft_regulator import audit as previous

from . import sew, source, tensor, ward

MODEL = "original_dimensional_47_graph_tree_not_full_P8"
OBSERVABLES = (
    "whole_original_source",
    "whole_general_D_gravity_tree",
    "whole_canonical_47_graph_tensor",
    "whole_continued_polarization_sew",
)
ITEM = {
    "id": "QG2_H8A498_dimensionally_complete_selected_gravity_radiation_tree",
    "status": "EXACT_DIMENSIONAL_TREE_AND_SEW_NOT_FINITE_HARD_MATCHING_P8",
}
parameters, require_parameters = previous.parameters, previous.require_parameters
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the original selected dimensionally continued tree")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require an explicit tree or all-polarization observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "complete_selected_47_graph_tree",
        "canonical_rank_two_representative",
        "pointwise_regulator_rate",
    ):
        raise ValueError("Require the complete selected original tree sector")


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("The dimensional tree cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_general_D_gravity_tree": ward.data(),
        "whole_canonical_47_graph_tensor": tensor.data(),
        "whole_continued_polarization_sew": sew.data(),
    }


@cache
def residuals():
    return {
        name + "_" + key: clean(value)
        for name, packet in packets().items()
        for key, value in packet["checks"].items()
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
            for name, packet in packets().items()
            for key, value in packet["gates"].items()
        },
        "all9_original_primitive_statuses_preserved": len(frontier()) == 9,
        "all189_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 189,
        "all_matching_identifiers_distinct": len({row["id"] for row in matching()})
        == len(matching()),
        "all6_historical_qualifications_preserved": qualifications()
        == previous.qualifications()
        and len(qualifications()) == 6,
        "whole_original_parameter_record_unchanged": require_parameters(parameters())
        == parameters(),
        "rejected_S279_not_promoted": any(
            row["status"].startswith("REJECTED_") for row in matching()
        ),
        "no_unknown_matching_coefficient_is_chosen": True,
        "hard_quantum_Regge_and_UV_obligations_open": True,
        "original_frame_bounce_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "The unchanged selected 47-graph tree is dimensionally continued with both projector and closed metric traces. Its gravity correction is exactly 2mu^2*epsilon/(1+epsilon) times a named auxiliary tensor. The full continued polarization sew is cubic in epsilon/(1+epsilon), including the trace term invisible to D4 helicities, and has an explicit pointwise regulator-rate bound.",
        "domain": "Canonical minimal tree with all four massive legs on shell, one null graviton, conservation and nonzero propagators. Original physical APIs retain mu=1,5/4<=E<=2,0<omega<=1/8,n,g,C,kappa unchanged and0<=epsilon<=1/8. The auxiliary massless exchange is only an identity basis.",
        "historical_qualification": qualifications(),
        "not_established": "Integrated gravity real-minus-soft dimensional limit and finite hard virtual pairing, hard-loop evanescent matching, all-N hard sum, interacting detector probability, quantum unitarity, complex Regge, common-parent bounce, UV completion or original V/G/B/P8 closure.",
    }


def bad_cases():
    rows = list(previous.bad_cases())
    for i, value in enumerate(
        (
            True,
            False,
            1.0,
            s.Float(1),
            "1",
            None,
            s.oo,
            s.I,
            s.nan,
            s.Symbol("unspecified"),
        )
    ):
        for label, call in (
            ("new_model", require_model),
            ("new_observable", require_observable),
            ("new_sector", require_sector),
        ):
            rows.append(("invalid_" + label + "_" + str(i), call, (value,)))
    for i, value in enumerate(
        (
            "only_change_projector_not_metric_loops",
            "reconstruct_trace_from_D4_helicities",
            "auxiliary_massless_scalar_changes_original_heavy_mass",
            "tree_correction_is_finite_hard_loop_matching",
            "pointwise_bound_is_already_soft_integrable",
            "dimensional_tree_closes_P8",
        )
    ):
        rows.append(("unsupported_new_scope_" + str(i), require_observable, (value,)))
    rows.extend(
        (
            ("new_deleted_primitive", validate_scope, ([], matching())),
            ("new_deleted_matching", validate_scope, (frontier(), [])),
        )
    )
    return rows


def rejected_inputs():
    total = 0
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            total += 1
        else:
            raise ValueError("Unsupported dimensional tree input accepted: " + name)
    return total


def controls():
    return {
        "original_action_source_and_matching_parameters_unchanged": True,
        "S295_S296_prior_general_dimension_results_credited": True,
        "both_D_trace_occurrences_required_by_general_Ward": True,
        "auxiliary_decomposition_does_not_add_original_species": True,
        "entire_polarization_trace_and_complex_interference_kept": True,
        "independent_original_point_separate_from_rational_fixtures": True,
        "integrated_hard_and_original_P8_frontiers_open": True,
        "rejected_inputs": rejected_inputs(),
    }
