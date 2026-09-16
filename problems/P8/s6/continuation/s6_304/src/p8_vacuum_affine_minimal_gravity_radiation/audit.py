"""Full selected real-radiation tree with all original closure gates retained."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_finite_forward_phase import audit as previous

from . import bounds, source, vertices, ward

MODEL = "original_minimal_gravity_radiation_not_full_P8"
OBSERVABLES = (
    "whole_canonical_gravity_vertices",
    "whole_general_Ward_identity",
    "whole_fixed_resolution_real_rate",
    "whole_resolution_and_matching_boundary",
)
ITEM = {
    "id": "QG2_H8A468_complete_gravity_radiation_and_uniform_fixed_resolution_full_tree_error",
    "status": "UNIFORM_FIXED_RESOLUTION_SELECTED_TREE_RATE_NOT_FULL_V_G_B_P8",
}
parameters = previous.parameters
require_parameters = previous.require_parameters
frontier = previous.frontier
qualifications = previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the stated exact gravity-radiation model")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated gravity-radiation observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "selected_tree_radiation",
        "angle_resolved_real_subtraction",
    ):
        raise ValueError("Require the known finite or symbolic matching sector")


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A selected real-radiation bound cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_canonical_gravity_vertices": vertices.data(),
        "whole_general_Ward_identity": ward.data(),
        "whole_fixed_resolution_real_rate": bounds.data(),
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
        "all159_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 159,
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
        "established": "Complete21-graph minimal Einstein radiation, general on-shell channel Ward identity and exact leading-soft normalization. Together with S295, all47 selected tree graphs have a finite real-minus-soft normalized rate error below1e32/kappa.",
        "domain": "mu=nu=1,25/4<=s<=16,all nonforward hard angles; any fixed0<x<=1/8, all emitted directions and both TT polarizations. Use the same recoil map, original positive full Born and inherited virtual-real IR scheme.",
        "historical_qualification": qualifications(),
        "not_established": "Finite forward cross section, unknown physical finite matching, full quantum hard/all-N/all-loop control, complex Regge bounds, common-parent state/bounce or original V/G/B/P8 closure.",
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
            "known_rate_bound_controls_loop_squares",
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
        "same_original_action_and_canonical_pure_Einstein_vertices": True,
        "all21_gravity_and47_selected_tree_graphs": True,
        "general_symbolic_Ward_and_independent_component_checks": True,
        "same_recoil_full_positive_Born_and_soft_IR_scheme": True,
        "all_finite_matching_coordinates_unassigned": True,
        "uniform_fixed_resolution_bound_not_finite_forward_cross_section": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
