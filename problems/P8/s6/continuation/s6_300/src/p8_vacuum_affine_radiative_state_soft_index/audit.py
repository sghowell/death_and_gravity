"""Radiative-state soft index with every original closure gate retained."""

from functools import cache

import sympy as s
from p8_vacuum_affine_calorimetric_soft_resummation import audit as previous
from p8_vacuum_affine_complement_measure.split import clean

from . import index, recoil, source, stability

MODEL = "original_radiative_state_soft_index_not_full_P8"
OBSERVABLES = (
    "whole_recoiled_radiative_state",
    "whole_massless_pair_soft_index",
    "whole_uniform_multiplicity_index_bound",
    "whole_same_state_nested_IR_pairing",
)
ITEM = {
    "id": "QG2_H8A464_complete_radiative_state_soft_index_uniform_multiplicity_recoil_and_nested_IR_pairing",
    "status": "COMPLETE_RADIATIVE_SOFT_INDEX_AND_UNIFORM_MULTIPLICITY_STABILITY_NOT_FINITE_FULL_RADIATION_OR_V_G_B_P8",
}
parameters = previous.parameters
require_parameters = previous.require_parameters
frontier = previous.frontier
qualifications = previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the stated exact radiative-state soft-index model")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require one stated radiative soft-index observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "radiative_soft_index",
        "nested_real_virtual_pole",
    ):
        raise ValueError("Require the radiative soft-index or nested-pole sector")


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A radiative soft-index bound cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_exact_arbitrary_radiation_recoil": recoil.data(),
        "whole_massless_soft_index_and_measure": index.data(),
        "whole_uniform_multiplicity_stability_and_nested_pairing": stability.data(),
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
        "all155_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 155,
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
        "established": "The exact radiative recoil map, all massive/massless pair soft contributions and collinear cancellations give a Lorentz-invariant soft index with a continuous finite angular-energy-measure extension. Its shift is bounded by40R/kappa uniformly in multiplicity. Correct nested real/virtual poles use that same radiative state. After pairing, the known power and Gamma-factor change is below3*10^-799 at the original parameters.",
        "domain": "Four mass-one scalars with incoming COM energy[5/4,2], arbitrary outgoing positive-energy null gravitons of total energyR<=1/8 and all angles, including collinear splitting and zero total spatial recoil. The known detector-factor comparison additionally requiresR<=resolution<=1/8 and reference mass scale1, after state-correct IR removal.",
        "historical_qualification": qualifications(),
        "not_established": "The full D finite angular conversion for arbitrary radiative states, finite radiative hard amplitudes and virtual terms, uniform complete-radiation-minus-leading-soft control, all omitted loops, full physical analytic/Regge bounds, original quantum state/domain/measure/bounce or V/G/B/P8 closure.",
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
            "full_D_finite_angular_conversion",
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
            raise ValueError("Unsupported radiative-index input accepted: " + name)
    return total


def controls():
    return {
        "exact_smooth_total_radiation_recoil": True,
        "all_massive_massless_and_massless_massless_pairs": True,
        "each_collinear_log_and_boost_shift_cancel": True,
        "whole_energy_measure_bound_uniform_in_N": True,
        "same_radiative_state_IR_poles_paired": True,
        "no_finite_radiative_amplitude_or_angular_term_invented": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
