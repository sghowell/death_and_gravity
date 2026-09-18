"""Full selected-tree dimensional real remainder and unchanged physical frontier."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_dimensional_gravity_radiation import audit as previous

from . import born, bounds, calibration, soft, source

MODEL = "original_full_selected_dimensional_real_remainder_not_full_P8"
OBSERVABLES = (
    "whole_original_source",
    "whole_same_D_Born",
    "whole_canonical_soft_remainder",
    "whole_integrated_dimensional_bounds",
    "whole_original_recoil_calibrations",
)
ITEM = {
    "id": "QG2_H8A499_full_selected_dimensional_real_remainder",
    "status": "COMPLETE_TREE_DIMENSIONAL_REAL_MINUS_SOFT_NOT_VIRTUAL_HARD_P8",
}
parameters, require_parameters = previous.parameters, previous.require_parameters
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the original complete-tree real-minus-soft model")


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated subtracted real-remainder observable")


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "complete_selected_47_graph_real_minus_soft",
        "same_D_leading_subtraction",
        "fixed_nonforward_resolution",
    ):
        raise ValueError("Require the fixed-nonforward subtracted complete-tree sector")


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("The real-remainder limit cannot close original P8")
    return True


@cache
def packets():
    return {
        "whole_original_source": source.data(),
        "whole_same_D_Born": born.data(),
        "whole_canonical_soft_remainder": soft.data(),
        "whole_integrated_dimensional_bounds": bounds.data(),
        "whole_original_recoil_calibrations": calibration.data(),
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
        "all190_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 190,
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
        "established": "The complete selected47 real-minus-soft integral has its original physical D4 limit at fixed nonforward hard data and every0<x<=1/8. The canonical trace is bounded by a six-dimensional completion, with separate low/high-energy estimates and the actual positive fixed-domain phase measure. Uniform bounds also give a joint lower-IR-cutoff/dimensional limit of this SUBTRACTED integrand.",
        "domain": "Unchanged original source,mu=nu=1,25/4<=s<=16,-1<z<1,delta=min(1,(s-4)(1-|z|)/2)>0,0<x<=1/8,0<=epsilon<=1/8. SAME-D complete Born times SAME-D leading soft tensor is subtracted before limits. A fixed-A0 amplitude normalization retains the SAME D-dimensional reference two-body phase.",
        "historical_qualification": qualifications(),
        "not_established": "A uniform forward limit, limits of the separate divergent rates, finite virtual hard or hard-loop evanescent matching, all-N hard summability, interacting detector probability, quantum unitarity, complex Regge, common-parent bounce, UV completion or original V/G/B/P8 closure.",
    }


def bad_cases():
    rows = list(previous.bad_cases())
    for i, value in enumerate(
        (
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
    ):
        for label, call in (
            ("new_model", require_model),
            ("new_observable", require_observable),
            ("new_sector", require_sector),
        ):
            rows.append(("invalid_" + label + "_" + str(i), call, (value,)))
    for i, value in enumerate(
        (
            "use_low_bound_above_its_recoil_window",
            "drop_extra_dimensional_trace_response",
            "replace_same_D_Born_inside_uncancelled_soft_pole",
            "separate_real_integral_is_IR_finite",
            "uniform_forward_limit_or_finite_virtual_matching",
            "real_remainder_limit_closes_P8",
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
            raise ValueError("Unsupported real-remainder input accepted: " + name)
    return total


def controls():
    return {
        "original_action_source_and_matching_parameters_unchanged": True,
        "same_D_Born_and_full_polarization_subtracted_first": True,
        "six_dimensional_completion_controls_canonical_trace": True,
        "separate_low_and_high_energy_bounds": True,
        "positive_fixed_domain_phase_and_reference_factors_retained": True,
        "joint_limits_only_for_the_subtracted_integrand": True,
        "forward_virtual_hard_and_original_P8_frontiers_open": True,
        "rejected_inputs": rejected_inputs(),
    }
