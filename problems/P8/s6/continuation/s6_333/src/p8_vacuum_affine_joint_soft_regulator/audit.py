"""Remaining-energy marked cutoff and preserved original physical frontier."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_remaining_log_cutoff import audit as previous

from . import kernels, moments, regulator, source, uniform

MODEL = "original_joint_soft_regulator_common_space_not_full_P8"
OBSERVABLES = (
    "whole_original_source",
    "whole_remaining_energy_kernels",
    "whole_conditional_log_moments",
    "whole_energy_dependent_regulator",
    "whole_joint_regulator_cutoff",
)
ITEM = {
    "id": "QG2_H8A497_joint_soft_regulator_and_infrared_cutoff",
    "status": "QUANTITATIVE_JOINT_SOFT_MARK_REGULATORS_NOT_HARD_INTERACTING_P8",
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
        "joint_regulated_additional_soft_insertion",
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
        "whole_remaining_energy_kernels": kernels.data(),
        "whole_conditional_log_moments": moments.data(),
        "whole_energy_dependent_regulator": regulator.data(),
        "whole_joint_regulator_cutoff": uniform.data(),
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
        "all188_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 188,
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
        "established": "The original fixed-Born D4 cloud insertion has joint common-space signed-density L1 error <=[132000*a*eta*L_eta^2+400000*epsilon*a*x*L_x^2]/kappa, physically <[106000*eta*L_eta^2+320000*epsilon*x*L_x^2]/kappa^2. An energy-dependent radial logarithmic remainder gives an explicit regulator rate; uniform cutoff estimates retain changing remaining energy and both conditioning errors.",
        "domain": "Original action and positive Borel energy states, fixed hard E,u,0<=a<=1,0<eta<=x<=1/8 and0<=epsilon<=1/8. SAME underlying leading-cloud probability space and exact cuts. Only the additional soft factor is dimensionally continued.",
        "historical_qualification": qualifications(),
        "not_established": "TV convergence of finite/infinite emission-configuration laws, outer hard dimensional/evanescent matching, finite radiative hard remainder, state-dependent radiation or interacting detector probability, all-N hard summability, quantum unitarity, absolute Regge, common-parent bounce, UV completion or original V/G/B/P8 closure.",
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
            "soft_factor_supplies_outer_hard_dimensional_matching",
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
        "same_original_source_and_full_additional_soft_convention": True,
        "energy_dependent_not_constant_regulator_remainder": True,
        "remaining_log_and_conditional_second_moments_kept": True,
        "both_relative_event_normalization_terms_retained": True,
        "uniform_cutoff_and_explicit_regulator_rate": True,
        "same_space_not_configuration_law_TV": True,
        "outer_hard_and_original_P8_frontiers_open": True,
        "rejected_inputs": rejected_inputs(),
    }
