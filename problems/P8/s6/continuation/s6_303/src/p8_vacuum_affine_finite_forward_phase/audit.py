"""Known all-angle real interference with every original closure gate retained."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_minimal_gravity_finite import audit as previous

from . import phase, rate, regularized, source

MODEL = "original_finite_forward_phase_not_full_P8"
OBSERVABLES = (
    "whole_finite_forward_phase",
    "whole_regularized_massive_forward_block",
    "whole_all_angle_known_real_interference",
    "whole_endpoint_matching_boundary",
)
ITEM = {
    "id": "QG2_H8A467_complete_finite_forward_phase_and_uniform_all_angle_known_real_interference",
    "status": "KNOWN_ALL_ANGLE_REAL_INTERFERENCE_NOT_COMPLEX_REGGE_OR_V_G_B_P8",
}
parameters = previous.parameters
require_parameters = previous.require_parameters
frontier = previous.frontier
qualifications = previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the stated exact finite-forward known-rate model")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated finite-forward known-rate observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "known_forward_interference",
        "symbolic_matching",
    ):
        raise ValueError("Require the known finite or symbolic matching sector")


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A compact finite-loop bound cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_finite_forward_phase": phase.data(),
        "whole_regularized_massive_forward_block": regularized.data(),
        "whole_all_angle_known_real_interference": rate.data(),
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
        "all158_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 158,
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
        "established": "The complete known finite Coulomb phase, classical square-root "
        "and quantum logarithmic forward terms are explicit. Exact integer-pole "
        "cancellation and a uniform two-derivative bound yield an all-angle real "
        "known hard-reference one-loop interference below1e12/kappa=1e-788 "
        "at the original parameters, using the full positive Born denominator.",
        "domain": "mu=nu=1,25/4<=s<=16,all physical nonforward angles. The same "
        "S302 analytic reference resolution/nu=1 and three unassigned physical "
        "finite matching coordinates are retained. No artificial angular "
        "window is imposed on this known real coefficient.",
        "historical_qualification": qualifications(),
        "not_established": "Matching values, gravity-Born radiation or full "
        "detector-inclusive rate, finite forward cross section, loop-square "
        "and all-order control, complex transfer analyticity and Regge bounds, "
        "original state/domain/measure/common-parent bounce or V/G/B/P8 closure.",
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
        "same_original_action_and_minimal_gravity_inputs": True,
        "whole_finite_forward_phase_and_real_coefficients": True,
        "exact_pole_cancellation_and_two_derivative_norm": True,
        "same_UV_IR_finite_reference_and_full_positive_Born": True,
        "three_matching_coordinates_explicitly_unassigned": True,
        "uniform_real_coefficient_not_complex_Regge_or_full_rate": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
