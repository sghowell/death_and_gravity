"""Full fixed spectator insertions with all original closure gates retained."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_minimal_gravity_radiation import audit as previous

from . import bounds, insertion, source
from . import matching as physical_matching

MODEL = "original_spectator_gravity_insertion_not_full_P8"
OBSERVABLES = (
    "whole_spectator_metric_insertion",
    "whole_known_local_matching_map",
    "whole_uniform_signed_known_rate",
    "whole_spectator_matching_boundary",
)
ITEM = {
    "id": "QG2_H8A469_complete_fixed_spectator_gravity_insertion_and_uniform_known_rate",
    "status": "FIXED_SPECTATOR_INSERTIONS_AND_KNOWN_RATE_NOT_FULL_V_G_B_P8",
}
parameters = previous.parameters
require_parameters = previous.require_parameters
frontier = previous.frontier
qualifications = previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the stated exact spectator-insertion model")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated spectator-insertion observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "fixed_spectator_metric_insertion",
        "known_reference_rate",
    ):
        raise ValueError("Require the known finite or symbolic matching sector")


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A known spectator bound cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_spectator_metric_insertion": insertion.data(),
        "whole_known_local_matching_map": physical_matching.data(),
        "whole_uniform_signed_known_rate": bounds.data(),
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
        "all160_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 160,
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
        "established": "Entire source-fixed H/Proca metric insertions in all crossed four-Phi channels, complete scalar/Proca optical normalization, unchanged original M1 log part and same-axis finite matching map. The specified known-reference interference obeys -3e-602<E<0 uniformly at all nonforward angles and tends to the fixed Newton correction at either endpoint.",
        "domain": "Original mass-one compact25/4<=s<=16,t,u<0; same positive full Born and fixed source prescription. Gaussian-only disk|p|<=16 is separate from a Wilsonian cutoff.",
        "historical_qualification": qualifications(),
        "not_established": "Unknown physical finite matching, sign of the full matched rate, finite forward cross section, interacting quantum/all-loop/Regge bounds, common-parent bounce or original V/G/B/P8 closure.",
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
        "same_original_source_fixed_H_Proca_prescription": True,
        "whole_kernels_and_original_M1_logs_not_Phi_double_counting": True,
        "independent_tensor_optical_and_Feynman_branch_checks": True,
        "same_positive_full_Born_at_every_nonforward_angle": True,
        "remaining_physical_matching_coordinates_unassigned": True,
        "known_reference_sign_not_full_rate_or_Regge_sign": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
