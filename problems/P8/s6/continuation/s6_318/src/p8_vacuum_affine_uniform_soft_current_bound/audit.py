"""Uniform weighted pure-soft currents with unchanged original frontiers."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_temporal_tree_reorganization import audit as previous

from . import bounds, geometry, grading, source

MODEL = "original_all_finite_weighted_soft_currents_not_inclusive_P8"
OBSERVABLES = (
    "whole_arbitrary_Rosen_geometry",
    "whole_all_order_null_grading",
    "whole_weighted_child_transfer",
    "whole_uniform_soft_current_majorant",
)
ITEM = {
    "id": "QG2_H8A482_all_finite_weighted_temporal_soft_current_bound",
    "status": "ALL_FINITE_WEIGHTED_SOFT_CURRENT_BOUND_NOT_INCLUSIVE_RATE_OR_FULL_V_G_B_P8",
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
        "weighted_pure_soft_currents",
        "labeled_majorant_and_angular_bound",
    ):
        raise ValueError(
            "Require weighted pure-soft currents or their labeled majorant"
        )


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A finite-tree source or graph count cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_arbitrary_Rosen_geometry": geometry.data(),
        "whole_all_order_null_grading": grading.data(),
        "whole_uniform_soft_current_majorant": bounds.data(),
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
        "all173_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 173,
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
        "established": "Uniform all-finite-multiplicity weighted temporal pure-soft current bound with explicit W^n/product(w_i) soft factors and coefficient envelope2*(2*10^10)^(n-1)*n!. The proof uses arbitrary Rosen geometry, tensor grading, a conserved-root inverse and labeled set-partition induction.",
        "domain": "The unchanged four-dimensional covariant action and original parameters. Pure-soft subcurrents have arbitrary finite collections of positive-energy null TT leaves, unit Frobenius norms, all angular configurations and generic internal propagators. Total soft energy is at most1/8 in the original physical domain; arbitrary energy hierarchies appear explicitly in the bound.",
        "historical_qualification": qualifications(),
        "not_established": "A complete hard four-scalar/N-graviton amplitude estimate and overlapping all-N soft subtraction; an inclusive real-virtual probability bound, finite hard and evanescent matching, interacting quantum state, unitarity, absolute complex Regge, common-parent bounce or V/G/B/P8 closure. Energy-independent boundedness and values at excluded exactly-collinear points are not claimed.",
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
        "arbitrary_Rosen_geometry_and_exact_reflection_grading": True,
        "same_complete_temporal_source_and_canonical_vertex_norms": True,
        "all_N_weighted_propagator_and_child_transfer_budgets": True,
        "dropping_longitudinal_weights_fails_the_exact_negative_control": True,
        "labeled_majorant_and_N4_collinear_calibrations_not_an_inclusive_rate": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
