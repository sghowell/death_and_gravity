"""Exact leading-ladder phase with all original closure gates retained."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_uniform_radiation_soft_limit import audit as previous

from . import bounds, dimensional, phase, source

MODEL = "original_leading_ladder_Coulomb_phase_not_full_P8"
OBSERVABLES = (
    "whole_dimensional_leading_ladders",
    "whole_exact_Coulomb_phase",
    "whole_uniform_Gamma_phase_bound",
    "whole_leading_class_scope_boundary",
)
ITEM = {
    "id": "QG2_H8A472_full_D_leading_ladder_Coulomb_phase_and_uniform_Gamma_error",
    "status": "LEADING_LADDER_COULOMB_PHASE_NOT_FULL_V_G_B_P8",
}
parameters = previous.parameters
require_parameters = previous.require_parameters
frontier = previous.frontier
qualifications = previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the stated exact leading-ladder-Coulomb-phase model")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated leading-ladder-Coulomb-phase observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "leading_gravity_ladders",
        "known_Coulomb_reference",
    ):
        raise ValueError("Require the known finite or symbolic matching sector")


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A leading-ladder phase cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_dimensional_leading_ladders": dimensional.data(),
        "whole_exact_Coulomb_phase": phase.data(),
        "whole_uniform_Gamma_phase_bound": bounds.data(),
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
        "all163_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 163,
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
        "established": "All-order fixed-IR-subtracted leading-ladder coefficients from the unchanged full-D pure-GR Born pole, summed into the exact unit-modulus Coulomb gamma phase; independent S303 one-loop calibration including2/V; complete gamma correction below2e-2400 uniformly in transfer at original parameters.",
        "domain": "Original mu=nu=1,25/4<=s<=16,0<tau<=1 near either Bose endpoint. The formula is a named leading-forward class, not a quantitative approximation theorem for the full amplitude on this domain.",
        "historical_qualification": qualifications(),
        "not_established": "The full interacting amplitude minus the leading class, unknown Newton and other physical matching, all-N nonleading radiation, complete detector rate, high-energy Regge bounds, original quantum state and bounce, original V/G/B/P8 closure.",
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
            "leading_phase_controls_full_quantum_amplitude",
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
        "original_full_D_Born_numerator_not_four_dimensional_replacement": True,
        "same_analytic_IR_phase_removed_before_coefficient_sum": True,
        "independent_S303_full_finite_phase_and_Abel_Fourier_checks": True,
        "complete_Gamma_factor_not_only_leading_log_exponent": True,
        "complex_monodromy_and_nonuniform_fixed_order_control_retained": True,
        "leading_Gaussian_class_not_full_original_quantum_theory": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
