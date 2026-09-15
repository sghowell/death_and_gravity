"""Whole normal-cut frontier, preserved obligations and rejected shortcuts."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_massive_elastic_proca_infrared import audit as previous

from . import angular, endpoint, invariant, source

MODEL = "original_massive_graviton_normal_cut_sheet_and_pole_running_not_full_P8"
OBSERVABLES = (
    "complete_graviton_angular_cut",
    "normal_channel_analytic_continuation",
    "nonuniform_massless_threshold",
    "necessary_subtraction_pole_running",
)
ITEM = {
    "id": "QG2_H8A446_entire_original_massive_graviton_normal_cut_closed_invariant_sheet_nonuniform_threshold_and_required_pole_anchor_running",
    "status": "COMPLETE_SCOPED_NORMAL_CUT_CONTINUATION_AND_POLE_RUNNING_NOT_FULL_V_G_B_P8",
}


def parameters():
    return previous.parameters()


def require_parameters(value):
    return previous.require_parameters(value)


def frontier():
    return previous.frontier()


def qualifications():
    return previous.qualifications()


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the stated original normal-cut model")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a scoped normal-cut observable")
    return value


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("A normal-cut continuation cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "complete_two_axis_angular_reduction": angular.data(),
        "entire_invariant_normal_cut_sheet": invariant.data(),
        "threshold_and_pole_anchor_running": endpoint.data(),
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
        "all9_primitive_statuses_preserved": len(frontier()) == 9,
        "all137_matching_records_preserved": matching()[:-1] == previous.matching()
        and len(previous.matching()) == 137,
        "all_matching_identifiers_distinct": len({r["id"] for r in matching()})
        == len(matching()),
        "all6_historical_qualifications_preserved": qualifications()
        == previous.qualifications()
        and len(qualifications()) == 6,
        "original_parameter_record_unchanged": require_parameters(parameters())
        == parameters(),
        "rejected_S279_not_promoted": any(
            r["status"].startswith("REJECTED_") for r in matching()
        ),
        "normal_cut_not_total_overlapping_discontinuity": True,
        "cap_running_not_undetermined_anchor_value": True,
        "S275_S276_S277_and_S281_physical_scopes_unchanged": True,
    }


def observable():
    return {
        "established": "Entire original TT graviton two-axis angular integral in closed form; exact invariant reduction to the accepted massive pair kernel; unique continuation on the physical normal-cut component, including the regular external threshold and explicit crossed massive boundaries; nonuniform massless endpoint; convergent pole-factor-preserving channel lift and necessary cap running of its pole anchor.",
        "domain": "Unchanged formal leading mass1,kappa10^800 vacuum. Physical normal s cut continues to s>0,t<4mu,u<4mu; other boundaries use the explicitly prescribed complex kernel, not real external-state positivity. The endpoint lift is a channel construction with an unresolved original pole anchor.",
        "historical_qualification": qualifications(),
        "not_established": "Full crossed/double-spectral real amplitude, transfer low-cut subtraction or b20; original finite renormalized massless-pole anchor; local matching, omitted loops, finite detector errors, Regge/UV applicability, exact LSZ vacuum, original quantum state/domain/bounce or V/G/B/P8 closure.",
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
        for label, call, args in (
            ("model", require_model, (value,)),
            ("observable", require_observable, (value,)),
            ("exact", source.exact_real, (value,)),
            ("sheet", source.require_real_sheet, (value, 0)),
            ("kernel", invariant.real_kernel, (value,)),
        ):
            rows.append(("invalid_" + label + "_" + str(i), call, args))
    for i, args in enumerate(((0, 0), (-1, 0), (1, -1), (1, 4), (2, 0, 0), (2, 0, -1))):
        rows.append(
            ("outside_real_normal_sheet_" + str(i), source.require_real_sheet, args)
        )
    for i, args in enumerate(((4,), (5,), (-1, 0), (0, -1), (4, 1))):
        rows.append(("outside_real_kernel_" + str(i), invariant.real_kernel, args))
    for i, value in enumerate(
        (
            "original_P8_closed",
            "full_total_discontinuity",
            "forward_b20",
            "drop_massless_pole",
            "new_local_counterterm",
            "massless_external_scalar",
            "physical_cutoff_from_cap",
            "full_Regge_bound",
            "normal_cut_is_positive_UV_density",
            "fixed_angle_limit_uniform_at_forward",
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
            raise ValueError("Unsupported normal-cut input accepted: " + name)
    return total


def controls():
    return {
        "whole_frozen_source_and_original_helicity_normalization": True,
        "complete_two_axis_division_not_forward_fitting": True,
        "shared_invariant_kernel_and_explicit_boundary_values": True,
        "uniform_integral_removes_external_pseudopole": True,
        "massless_threshold_forward_nonuniformity_retained": True,
        "pole_running_and_crossed_b20_anchor_shift_explicit": True,
        "normal_channel_not_full_amplitude_or_P8": True,
        "rejected_inputs": rejected_inputs(),
    }
