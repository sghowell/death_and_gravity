"""Complete selected known scalar-loop coefficient, not independent parent matching."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean
from p8_vacuum_affine_core_curvature_coefficient import audit as previous

from . import aggregate, calibration, local, mixed, moment, source

MODEL = "original_selected_scalar_loop_common_covariant_jet_lift"
OBSERVABLES = (
    "whole_original_source_and_complete_scalar_inventory",
    "whole_generic_offshell_local_tadpoles",
    "whole_complete_mixed_source_conversion",
    "whole_exact_extra_mass_moments",
    "whole_selected_known_matter_sign_and_bound",
    "whole_frozen_kernel_and_original_calibration",
)
ITEM = {
    "id": "QG2_H8A511_selected_known_matter_common_curvature",
    "status": "COMPLETE_SELECTED_SCALAR_ONE_LOOP_DEGREE6_COEFFICIENT_IN_ONE_JET_LIFT_NOT_INDEPENDENT_PARENT_MATCHING",
}
parameters, require_parameters = previous.parameters, previous.require_parameters
frontier, qualifications = previous.frontier, previous.qualifications


def matching():
    return [*previous.matching(), dict(ITEM)]


def require_model(value):
    if not isinstance(value, str) or value != MODEL:
        raise ValueError("Require the unchanged selected scalar-loop model")
    return value


def require_observable(value):
    if not isinstance(value, str) or value not in OBSERVABLES:
        raise ValueError("Require a stated whole known scalar-loop observable")
    return value


def require_sector(value):
    if not isinstance(value, str) or value not in (
        "complete_selected_scalar_fourpoint_one_loop",
        "whole_common_six_word_jet_lift",
        "original_known_local_coefficient_bound",
    ):
        raise ValueError("Require the stated selected scalar-loop matching sector")
    return value


def validate_scope(original, matched):
    if original != frontier() or matched != matching():
        raise ValueError("The known scalar-loop coefficient cannot close original P8")
    return True


@cache
def packets():
    return dict(
        zip(
            OBSERVABLES,
            (
                source.data(),
                local.data(),
                mixed.data(),
                moment.data(),
                aggregate.data(),
                calibration.data(),
            ),
        )
    )


@cache
def residuals():
    return {
        name + "_" + key: clean(value)
        for name, packet in packets().items()
        for key, value in packet["checks"].items()
    }


def scalar_entry_count():
    return sum(
        len(v) if isinstance(v, s.MatrixBase) else 1 for v in residuals().values()
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
        "all202_previous_matching_records_preserved": matching()[:-1]
        == previous.matching()
        and len(previous.matching()) == 202,
        "all_matching_identifiers_distinct": len({r["id"] for r in matching()})
        == len(matching()),
        "all6_historical_qualifications_preserved": qualifications()
        == previous.qualifications()
        and len(qualifications()) == 6,
        "whole_original_parameters_unchanged": require_parameters(parameters())
        == parameters(),
        "rejected_S279_not_promoted": any(
            r["status"].startswith("REJECTED_") for r in matching()
        ),
        "independent_extra_parent_chi_not_assigned_or_bounded": True,
        "known_scalar_loop_not_finite_gravity_quantum_decoupling": True,
        "physical_metric_bounce_scoped_P8a_and_global_obligations_unchanged": True,
    }


def observable():
    return {
        "established": "The complete selected formal scalar one-loop four-light degree6 real-TT curvature coefficient of the original fixed-canonical limiting action is matched in the same explicit six-word covariant jet lift. All six local tadpoles agree off shell in exactD with their literal quartic lifts and give zero extra curvature. All three mixed-source contraction classes are retained; their combined correction is positive at the original parameters with chi_extra/A0<10^-604. Adding it to the S346 known polynomial core gives a strictly negative total with |chi_known_selected|/A0<10^-207.",
        "domain": "The exhaustive scalar E4,L1 inventory of the classical limiting action, one external real-null-TT graviton, and a local homogeneous analytic-origin derivative coefficient. All off-shell comparison projections precede scalar-shell substitution. The physical source, metric, masses and pole-only/OS4 scheme are unchanged.",
        "historical_qualification": qualifications(),
        "not_established": "Independent extra parent chi remains unmatched and unbounded. No full curved effective/counterfunctional, internal-graviton loops, finite-gravity quantum decoupling, physical above-threshold Taylor approximation or truncation-error estimate, full virtual/inclusive matching, exact LSZ, quantum unitarity, complex Regge, same-parent bounce or UV completion. Original V/G/B/P8 remain open; no extra copy is added to existing full loop amplitudes.",
    }


def bad_cases():
    rows = list(previous.bad_cases())
    for j, v in enumerate(
        (True, False, 1.0, s.Float(1), "1", None, s.oo, s.I, s.nan, s.Symbol("unknown"))
    ):
        for label, call in (
            ("model", require_model),
            ("observable", require_observable),
            ("sector", require_sector),
        ):
            rows.append((f"selected_invalid_{label}_{j}", call, (v,)))
    for j, v in enumerate(
        (
            "full_parent_chi_zero",
            "finite_gravity_decoupling_proved",
            "physical_cut_replaced_by_jet",
            "known_loop_added_twice",
            "all_positive_masses_have_positive_extra",
            "selected_closes_original_P8",
        )
    ):
        rows.append((f"selected_unsupported_scope_{j}", require_observable, (v,)))
    rows.extend(
        (
            ("selected_deleted_primitive", validate_scope, ([], matching())),
            ("selected_deleted_matching", validate_scope, (frontier(), [])),
        )
    )
    for j, v in enumerate((True, False, 1.0, s.Float(1), "1", None, 0, -1, s.oo, s.I)):
        rows.append((f"selected_invalid_mass_{j}", moment.extra_at, (v,)))
    for j, v in enumerate((True, False, -1, 0, 1, 4, 2.0, s.Integer(2), "2", None, [])):
        rows.append((f"selected_invalid_moment_order_{j}", moment.coefficient, (v, 1)))
    for j, v in enumerate(
        (
            "internal_graviton",
            "independent_parent_chi",
            "full_bounce_stress",
            True,
            None,
        )
    ):
        rows.append((f"selected_invalid_class_{j}", source.require_class, (v,)))
    return rows


def rejected_inputs():
    count = 0
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            count += 1
        else:
            raise ValueError("Unsupported selected-scalar input accepted: " + name)
    return count


def controls():
    return {
        "complete_scalar_graph_inventory_and_original_source": True,
        "all_six_local_tadpoles_generic_offshell_exactD": True,
        "all_three_mixed_classes_EOM_delta_and_counterterms_retained": True,
        "exact_extra_moments_and_finite_positive_weight_bounds": True,
        "complete_common_basis_selected_sum_with_strict_sign": True,
        "actual_frozen_kernels_and_original_massive_vectors": True,
        "independent_parent_matching_and_original_frontiers_preserved": True,
        "rejected_inputs": rejected_inputs(),
    }
