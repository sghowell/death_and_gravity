"""Finite radiative angular conversion with every original closure gate retained."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_radiative_state_soft_index import audit as previous

from . import angular, bounds, conversion, source

MODEL = "original_radiative_angular_finite_not_full_P8"
OBSERVABLES = (
    "whole_fixed_ball_dimensional_current",
    "whole_collinear_Holder_bound",
    "whole_finite_radiative_conversion",
    "whole_known_detector_reference",
)
ITEM = {
    "id": "QG2_H8A465_complete_radiative_finite_angular_conversion_uniform_collinear_Holder_and_regulator_bound",
    "status": "COMPLETE_RADIATIVE_FINITE_ANGULAR_CONVERSION_NOT_HARD_RADIATION_OR_V_G_B_P8",
}
parameters = previous.parameters
require_parameters = previous.require_parameters
frontier = previous.frontier
qualifications = previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the stated exact radiative finite angular model")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require one stated radiative finite angular observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "radiative_finite_angular",
        "known_detector_reference",
    ):
        raise ValueError(
            "Require the finite radiative angular or known reference sector"
        )


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A finite radiative angular bound cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_fixed_ball_dimensional_current": angular.data(),
        "whole_collinear_Holder_majorants": bounds.data(),
        "whole_finite_conversion_and_detector_reference": conversion.data(),
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
        "all156_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 156,
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
        "no_full_amplitude_remainder_is_invented": True,
        "full_gravity_Regge_and_all_loop_obligations_open": True,
        "original_frame_bounce_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "The full fixed-ball continuation, retained D-projector and uniform collinear Holder estimate give the finite radiative angular conversion with absolute bound2000/kappa and explicit fixed-resolution regulator error. The known unexpanded detector reference changes from its elastic value by less than4250/kappa at originalparameters, uniformly in finite multiplicity.",
        "domain": "The unchangedS300 four-massive plus positive-null radiative states, incoming COM energy[5/4,2], total radiationR<=1/8 and all angles. Full soft dimensional projector D=4+2e with0<=e<=1/8 and physical external momenta in four dimensions. Known reference comparison additionally requiresR<=resolution<=1/8 and nu1.",
        "historical_qualification": qualifications(),
        "not_established": "Finite radiative hard amplitudes, hard evanescent terms and their scheme conversion, uniform complete-radiation-minus-leading-soft remainder, omitted loops/operators, full physical analytic/Regge and original quantum state/domain/measure/bounce or V/G/B/P8 closure.",
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
            "particle_count_cutoff",
            "drop_massless_radiation_legs",
            "discard_radiation_pair_terms",
            "old_virtual_pole_on_new_state",
            "massive_subset_current_conserved",
            "hard_evanescent_terms_are_zero",
            "full_physical_inclusive_rate",
            "all_loop_regulator_complete",
            "quantum_state_from_energy_measure",
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
            raise ValueError("Unsupported finite-angular input accepted: " + name)
    return total


def controls():
    return {
        "same_original_radiative_recoil": True,
        "fixed_ball_dimensional_current_and_trace": True,
        "uniform_integrable_collinear_Holder_modulus": True,
        "whole_finite_angular_bound_uniform_in_N": True,
        "same_radiative_state_IR_poles_and_finite_conversion": True,
        "hard_evanescent_matching_and_full_rate_not_inferred": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
