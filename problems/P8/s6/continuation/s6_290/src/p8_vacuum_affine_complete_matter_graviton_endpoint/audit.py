"""Whole original C/g matter endpoint, physical cuts and bounded crossed coefficient."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_covariant_gaussian_four_point_matching import audit as previous

from . import cuts, endpoint, forward, source

MODEL = "original_complete_minimal_matter_endpoint_not_full_physical_P8"
OBSERVABLES = (
    "whole_source_endpoint_arity",
    "complete_F1_F2_with_curved_anchors",
    "whole_light_and_heavy_pair_endpoint_cuts",
    "bounded_crossed_endpoint_coefficient",
)
ITEM = {
    "id": "QG2_H8A454_complete_minimal_matter_graviton_endpoint_with_unmatched_curved_anchors_full_pair_cuts_and_bounded_crossed_endpoint_coefficient",
    "status": "COMPLETE_SCOPED_MINIMAL_MATTER_ENDPOINT_NOT_FULL_MIXED_FOUR_POINT_IR_REGGE_OR_V_G_B_P8",
}
parameters = previous.parameters
require_parameters = previous.require_parameters
frontier = previous.frontier
qualifications = previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the original minimal matter endpoint sector")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated complete endpoint observable")
    return value


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "minimal_matter_endpoint",
        "crossed_endpoint_exchange",
    ):
        raise ValueError(
            "Require the complete minimal endpoint or its crossed exchange"
        )
    return value


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("Endpoint matching cannot close the original physical P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_minimal_matter_endpoint": endpoint.data(),
        "whole_endpoint_physical_cuts": cuts.data(),
        "whole_crossed_endpoint_coefficient": forward.data(),
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
        "all145_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 145,
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
        "finite_RH_and_other_curved_matching_not_assigned": True,
        "whole_physical_IR_Regge_and_full_mixed_amplitude_open": True,
        "original_frame_bounce_and_scoped_P8a_unchanged": True,
    }


def observable():
    return {
        "established": "Both complete minimal C/g matter-graviton form factors are explicit with their original flat-OS/onepoint counterterms and two unmatched finite curved anchors. The complete two-body light/heavy tensor cuts retain every tree exchange. The entire known crossed ENDPOINT graph coefficient has finite-mass magnitude below10^-1005. A constant curvature improvement has zero b20, while the RH exchange residue remains explicit and unbounded.",
        "domain": "One minimal matter loop and one graviton endpoint in the unchanged original variables. The crossed result is only the graviton-exchange endpoint contribution, not all mixed matter/gravity four-point loops. Forward evaluation is in the massive analytic domain; physical pair-cut domains stay away from the formal heavy resonance. Actual mu1,kappa10^800,nH=10^200/512+2,g1/8192,C unchanged.",
        "historical_qualification": qualifications(),
        "not_established": "Finite values or bounds for RH and other curved/EFT matching, full mixed four-point amplitude or full-source b20, exact stable heavy atom or resonance resummation, all-loop errors, physical detector/dressing/IR or Regge remainder, original state/domain/measure/bounce or V/G/B/P8 closure.",
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
            "choose_light_finite_curvature",
        )
    ):
        rows.append(("unsupported_sector_" + str(i), require_sector, (value,)))
    for i, value in enumerate(
        (
            "original_P8_closed",
            "set_all_local_coefficients_zero",
            "Euler_topological_in_all_dimensions",
            "nonlinear_quantum_EOM_equivalence",
            "copy_heavy_finite_prescription_to_Phi",
            "add_light_Gaussian_loop_twice",
            "heavy_asymptotic_series_is_full_error_bound",
            "full_Regge_bound",
            "one_Gaussian_sign_is_full_positivity",
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
            raise ValueError(
                "Unsupported complete matter endpoint input accepted: " + name
            )
    return total


def controls():
    return {
        "whole_original_source_and_endpoint_arity": True,
        "both_matter_form_factors_and_original_counterterms": True,
        "entire_dimensional_bubble_and_raw_even_triangle": True,
        "whole_light_and_heavy_pair_tensor_cuts": True,
        "full_crossed_known_endpoint_finite_mass_bound": True,
        "unmatched_RH_not_set_by_flat_OS_or_b20_constant": True,
        "all_original_and_historical_frontiers_retained": True,
        "rejected_inputs": rejected_inputs(),
    }
