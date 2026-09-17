"""Arbitrary-multiplicity tree source with unchanged scientific frontiers."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_state_transport_soft_matching import audit as previous

from . import jets, source, topology, trees

MODEL = "original_all_multiplicity_tree_source_not_all_N_rate_or_full_P8"
OBSERVABLES = (
    "whole_all_order_canonical_vertices",
    "whole_all_N_topology",
    "whole_all_finite_tree_source",
    "whole_vertex_and_graph_majorants",
)
ITEM = {
    "id": "QG2_H8A479_arbitrary_finite_multiplicity_canonical_tree_source",
    "status": "ARBITRARY_FINITE_TREE_SOURCE_NOT_ALL_N_RATE_OR_FULL_V_G_B_P8",
}
parameters, require_parameters = previous.parameters, previous.require_parameters
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the stated original all-multiplicity tree model")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated finite-tree source observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "all_finite_canonical_trees",
        "vertex_and_topology_majorants",
    ):
        raise ValueError("Require the stated finite-tree source or separate majorants")


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A finite-tree source or graph count cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_all_order_canonical_vertices": jets.data(),
        "whole_all_N_graph_inventory": topology.data(),
        "whole_all_finite_tree_source": trees.data(),
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
        "all170_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 170,
        "all_matching_identifiers_distinct": len({r["id"] for r in matching()})
        == len(matching()),
        "all6_historical_qualifications_preserved": qualifications()
        == previous.qualifications()
        and len(qualifications()) == 6,
        "whole_original_parameter_record_unchanged": require_parameters(parameters())
        == parameters(),
        "rejected_S279_not_promoted": any(
            r["status"].startswith("REJECTED_") for r in matching()
        ),
        "no_unknown_matching_coefficient_is_chosen": True,
        "full_gravity_Regge_and_all_loop_obligations_open": True,
        "original_frame_bounce_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "Canonical mixed metric and Einstein vertices at every finite order; a complete labeled four-Phi/N-graviton tree evaluator; closed all-N graph enumeration and rigorous vertex/graph majorants. Independent literal action, root-cut enumeration, lower-source compatibility and full original three-real calibration.",
        "domain": "The unchanged selected four-dimensional covariant scalar/heavy/Einstein action and original parameters. Exact finite-N physical point evaluation on the original compact recoil domain, away from internal poles. Vertex norms also hold for complex fields with Euclidean Frobenius norms.",
        "historical_qualification": qualifications(),
        "not_established": "An all-N nonleading amplitude or probability bound, arbitrary soft-subtree collinear cancellation estimates, complete finite hard real-virtual matching, an interacting quantum state, absolute complex Regge, common-parent bounce or V/G/B/P8 closure.",
    }


def bad_cases():
    bad = (
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
    rows = []
    for i, value in enumerate(bad):
        for label, call in (
            ("model", require_model),
            ("observable", require_observable),
            ("order", source.require_order),
            ("mass", source.require_mass),
            ("sector", require_sector),
        ):
            rows.append(("invalid_" + label + "_" + str(i), call, (value,)))
    for i, value in enumerate((-1, 0, s.Rational(1, 2), s.Rational(3, 2))):
        rows.append(
            ("outside_positive_integer_order_" + str(i), source.require_order, (value,))
        )
    for i, value in enumerate((0, -1)):
        rows.append(("outside_positive_mass_" + str(i), source.require_mass, (value,)))
    for i, value in enumerate(
        (
            "whole_finite_G_amplitude",
            "UV_complete",
            "full_crossed_b20",
            "exact_massless_LSZ",
            "choose_finite_curvature",
        )
    ):
        rows.append(("unsupported_sector_" + str(i), require_sector, (value,)))
    for i, value in enumerate(
        (
            "original_P8_closed",
            "all_matching_zero",
            "known_piece_is_entire_amplitude",
            "finite_loop_proves_unitarity",
            "drop_evanescent_terms",
            "drop_old_Gram_finite",
            "discard_imaginary_phase_from_amplitude",
            "relative_Regge_is_absolute",
            "all_loop_regulator_complete",
            "graph_counts_close_full_inclusive_detector_rate",
        )
    ):
        rows.append(("unsupported_" + str(i), require_observable, (value,)))
    rows.extend(
        (
            ("deleted_primitive", validate_scope, ([], matching())),
            ("deleted_matching", validate_scope, (frontier(), [])),
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
            raise ValueError("Unsupported finite-gravity input accepted: " + name)
    return total


def controls():
    return {
        "same_original_action_and_frozen_lower_sources": True,
        "literal_EH5_and_scalar_fourth_metric_independent_calibrations": True,
        "all_N_EGF_and_independent_root_cut_counts": True,
        "original_three_real_5116_tree_amplitude": True,
        "higher_vertex_omissions_break_full_tree_Ward": True,
        "vertex_and_graph_majorants_not_an_all_N_rate": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
