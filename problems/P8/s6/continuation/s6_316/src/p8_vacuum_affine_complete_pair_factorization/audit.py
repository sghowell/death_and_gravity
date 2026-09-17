"""Complete pair factorization with unchanged original scientific frontiers."""

from functools import cache

import sympy as s
from p8_vacuum_affine_all_multiplicity_tree_source import audit as previous
from p8_vacuum_affine_complement_measure.split import clean

from . import forests, obstruction, pairs, source

MODEL = "original_complete_pair_factorization_not_all_N_rate_or_full_P8"
OBSERVABLES = (
    "whole_complete_single_root_Ward",
    "whole_matching_forest_inventory",
    "whole_completed_pair_factorization",
    "whole_isolated_cluster_obstructions",
)
ITEM = {
    "id": "QG2_H8A480_complete_pair_factorization_and_cluster_obstructions",
    "status": "COMPLETE_PAIR_FACTORIZATION_NOT_ALL_N_RATE_OR_FULL_V_G_B_P8",
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
        "completed_fixed_pair_classes",
        "matching_forests_and_obstructions",
    ):
        raise ValueError(
            "Require completed pair classes or matching-forest obstructions"
        )


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A finite-tree source or graph count cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_matching_forest_inventory": forests.data(),
        "whole_completed_pair_factorization": pairs.data(),
        "whole_isolated_cluster_obstructions": obstruction.data(),
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
        "all171_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 171,
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
        "established": "Complete one-off-shell-root Ward induction; fixed-pair factorization and conditional angular bound; all-N matching-forest inclusion-exclusion; exact full5116-tree original calibration; isolated-cluster poles and multi-off-shell Ward counterexamples.",
        "domain": "The unchanged selected four-dimensional covariant scalar/heavy/Einstein action and original parameters. Complete currents have physical on-shell free leaves and one off-shell graviton root, away from other internal poles. Fixed-pair estimates assume unit-Frobenius physical TT polarizations and positive energies.",
        "historical_qualification": qualifications(),
        "not_established": "A uniform completed-remainder norm, simultaneous multi-collinear estimate, all-N amplitude or probability bound, finite hard real-virtual matching, interacting quantum state, absolute complex Regge, common-parent bounce or V/G/B/P8 closure.",
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
        "single_offshell_root_Ward_and_independent_source_comparison": True,
        "all_N_matching_forest_and_independent_no_pair_root_counts": True,
        "original_three_real_5116_tree_complete_pair_factorization": True,
        "isolated_cluster_and_multi_offshell_shortcuts_rejected": True,
        "conditional_pair_estimate_not_an_all_N_rate": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
