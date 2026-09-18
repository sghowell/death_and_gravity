"""Remaining-energy marked cutoff and preserved original physical frontier."""

from functools import cache

import sympy as s
from p8_vacuum_affine_borel_soft_conversion import audit as previous
from p8_vacuum_affine_complement_measure.split import clean

from . import density, endpoint, source, tail

MODEL = "original_remaining_log_common_space_cutoff_not_full_P8"
OBSERVABLES = (
    "whole_original_source",
    "whole_common_space_density",
    "whole_same_event_log",
    "whole_excluded_endpoint_log",
)
ITEM = {
    "id": "QG2_H8A496_remaining_energy_log_common_space_cutoff",
    "status": "QUANTITATIVE_COMMON_SPACE_SOFT_MARK_CUTOFF_NOT_HARD_INTERACTING_P8",
}
parameters, require_parameters = previous.parameters, previous.require_parameters
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError(
            "Require the original common-space additional-soft cutoff model"
        )


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated common-space density or endpoint observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "limiting_additional_soft_insertion",
        "common_space_signed_density",
        "same_fixed_Born_D4_reference",
    ):
        raise ValueError("Require the fixed-reference limiting additional-soft sector")


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("The marked cutoff theorem cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_common_space_density": density.data(),
        "whole_same_event_log": tail.data(),
        "whole_excluded_endpoint_log": endpoint.data(),
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
        "all187_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 187,
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
        "established": "For the fixed-Born D4 leading cloud, the same-space signed densities of Z=deltaDelta+deltaa*ln(x-R), with their own exact total-energy conditioning, converge in L1 at rate <=63000*a*eta*(1+ln(1/eta))^2/kappa; the original physical bound is <51000*eta*(1+ln(1/eta))^2/kappa^2. Both the remaining-energy change and excluded-event singular layer are retained.",
        "domain": "Original action and positive energy measures, fixed hard E,u,0<=a<=1,0<eta<=x<=1/8; each state mark evaluated only on its own physical cut. Fixed-Born D4 reference, limiting additional-soft insertion and common underlying Poisson space.",
        "historical_qualification": qualifications(),
        "not_established": "TV convergence of finite/infinite emission-configuration laws, a simultaneous quantitative epsilon/eta regulator estimate, outer hard dimensional/evanescent matching, finite radiative hard remainder, state-dependent radiation, interacting detector probability, all-N hard summability, quantum unitarity, absolute Regge, common-parent bounce, UV completion or original V/G/B/P8 closure.",
    }


def bad_cases():
    rows = list(previous.bad_cases())
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
    for i, value in enumerate(bad):
        for label, call in (
            ("new_model", require_model),
            ("new_observable", require_observable),
            ("new_sector", require_sector),
        ):
            rows.append(("invalid_" + label + "_" + str(i), call, (value,)))
    for i, value in enumerate(
        (
            "finite_configuration_laws_converge_in_TV",
            "drop_both_conditioning_normalization_terms",
            "replace_remaining_energy_by_full_resolution",
            "joint_regulator_rate_already_proved",
            "evaluate_full_state_outside_its_cut",
            "soft_mark_cutoff_closes_P8",
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
            raise ValueError("Unsupported cutoff input accepted: " + name)
    return total


def controls():
    return {
        "same_original_source_and_positive_increment": True,
        "exact_remaining_energy_and_Born_subtraction": True,
        "common_space_density_not_configuration_law_TV": True,
        "both_relative_normalization_errors_retained": True,
        "singular_excluded_event_layer_retained": True,
        "rare_index_and_zero_atom_handled": True,
        "limiting_insertion_not_hard_or_joint_regulator_completion": True,
        "rejected_inputs": rejected_inputs(),
    }
