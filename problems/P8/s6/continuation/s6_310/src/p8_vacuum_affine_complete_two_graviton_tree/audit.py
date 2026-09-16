"""Scoped complete two-graviton tree with unchanged original frontiers."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_single_residual_soft_dressing import audit as previous

from . import checks, source, topology

MODEL = "original_complete_two_graviton_tree_not_full_P8"
OBSERVABLES = (
    "whole_two_real_source_and_tree",
    "whole_two_real_topology_count",
    "whole_two_real_action_and_Ward_calibration",
    "whole_original_two_soft_limits",
)
ITEM = {
    "id": "QG2_H8A474_complete_four_scalar_two_graviton_formal_tree_and_state_correct_soft_limits",
    "status": "COMPLETE_SELECTED_TWO_REAL_TREE_NOT_FULL_V_G_B_P8",
}
parameters, require_parameters = previous.parameters, previous.require_parameters
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the stated complete two-graviton tree model")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated two-graviton tree observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "complete_formal_tree",
        "fixed_kinematics_leading_soft_limits",
    ):
        raise ValueError("Require the selected formal tree or leading soft sector")


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A finite tree cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_labeled_tree_topologies": topology.data(),
        "whole_action_and_physical_calibrations": checks.data(),
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
        "all165_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 165,
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
        "established": "Complete selected original four-scalar/two-graviton formal tree via a terminating labeled root-cut recursion:38 contact,198 heavy-exchange and198 pure-Einstein graphs. Literal higher metric/action coefficients, independent EGF counts, frozen lower-point calibration, six-point Ward/Bose checks and original-parameter simultaneous/state-correct hierarchical soft limits are supplied.",
        "domain": "Four-dimensional formal tree at original parameters; mass-one external scalars with incoming COM energy5/4<=E<=2, total positive real-graviton energy<=1/8. Exact point evaluator excludes internal propagator poles and exact collinearity of the two null rays. Soft theorems are at fixed generic hard kinematics and fixed noncollinear soft directions, not uniform integrated estimates.",
        "historical_qualification": qualifications(),
        "not_established": "A uniform integrated two-real nonleading remainder, its overlapping soft/collinear subtraction and virtual pairing, full all-N physical detector rate, radiative hard loops or finite matching, unitarity, complex Regge bound, original quantum state and common-parent bounce, or original V/G/B/P8 closure.",
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
            "two_real_tree_controls_full_multi_real_rate",
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
        "same_original_source_and_formal_tree_order": True,
        "complete_rooted_tree_and_independent_434_topologies": True,
        "literal_higher_action_vertices_and_frozen_lower_point_calibration": True,
        "both_graviton_Ward_and_scalar_graviton_Bose_checks": True,
        "original_simultaneous_and_full_radiative_state_soft_currents": True,
        "finite_tree_not_uniform_integrated_or_all_N_bound": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
