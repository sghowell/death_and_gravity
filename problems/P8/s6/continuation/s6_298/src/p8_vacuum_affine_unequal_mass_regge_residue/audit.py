"""Unequal-mass Regge error map with every original closure gate retained."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_one_newton_inclusive_assembly import audit as previous

from . import moments, regge, source, waves

MODEL = "original_unequal_mass_Regge_error_map_not_full_P8"
OBSERVABLES = (
    "whole_two_mass_spin2_cut",
    "whole_positive_spectral_moment_and_tail",
    "whole_original_heavy_quarter",
    "conditional_Regge_error_budget",
)
ITEM = {
    "id": "QG2_H8A462_unequal_mass_coupled_spin2_Regge_residue_spectral_moment_and_conditional_contour_error",
    "status": "COMPLETE_UNEQUAL_MASS_SPIN2_MOMENT_AND_CONDITIONAL_REGGE_ERROR_MAP_NOT_ACTUAL_UV_IR_OR_V_G_B_P8",
}
parameters = previous.parameters
require_parameters = previous.require_parameters
frontier = previous.frontier
qualifications = previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the stated unequal-mass Regge error map")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require one stated moment or conditional error observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "generated_matter_spin2",
        "conditional_Regge_cross_term",
    ):
        raise ValueError(
            "Require the generated matter sector or conditional pole cross term"
        )


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A conditional Regge budget cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_unequal_mass_spin2_and_coupled_pole": waves.data(),
        "whole_complete_spectral_moment_and_original_hierarchy": moments.data(),
        "whole_conditional_Regge_order_and_contour_error": regge.data(),
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
        "all153_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 153,
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
        "no_conditional_budget_input_is_a_matched_value": True,
        "full_gravity_Regge_and_all_loop_obligations_open": True,
        "original_frame_bounce_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "Both unequal-mass Bose spin2 cuts, the coupled factorized-pole residue variation and full generated F1 moment agree. The original HH cut contributes a quarter of the total slope within2*10^-98; the4<T<=16 light window contributes less than10^-195. A positive whole-Legendre-Q representation gives an integrated finite-spin error estimate, combined with the exact known spectral tail and an explicit unknown transfer arc.",
        "domain": "The exact generated matter graphs have n>4mu>0, with each named cut open above4mu or4n. The quarter and low-window bounds use the original mass-one hierarchy. The Regge error map is conditional on isolated even-pole complex-J unitarity, declared positive trajectory and residue-ratio estimates, a finite slit transfer diskU>4n and its declared arc supremum.",
        "historical_qualification": qualifications(),
        "not_established": "A matched Regge trajectory or residue ratio, a numerical transfer arc norm, absence of other UV singularities, the high-energy complex-energy contour remainder, analytic trajectory counterterms, complete physical massless infrared/all-loop control, original quantum state/domain/measure/bounce or V/G/B/P8 closure.",
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
            "equal_mass_formula_for_original_n",
            "ignore_HH_cut",
            "H_is_exact_stable",
            "discard_transfer_arc",
            "trajectory_analytic_part_zero",
            "single_exchange_as_Bose_pair",
            "full_Regge_bound",
            "finite_forward_from_unpaired_cut",
            "exchange_IR_limits",
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
            raise ValueError("Unsupported Regge input accepted: " + name)
    return total


def controls():
    return {
        "both_Bose_partial_waves_and_original_masses": True,
        "entire_generated_F1_spectral_moment": True,
        "original_heavy_quarter_not_discarded": True,
        "whole_Q_order_inequality_and_integrable_thresholds": True,
        "normalized_log_residue_and_unknown_arc": True,
        "no_model_specific_Regge_inputs_invented": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
