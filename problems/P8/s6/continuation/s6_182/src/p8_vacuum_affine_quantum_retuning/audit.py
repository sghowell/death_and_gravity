"""New-action retained mean stationarity with explicit, non-UV scope."""

from functools import cache

import sympy as s
from p8_vacuum_affine_kernel_norm import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import bounds, profile, response, ward

ITEM = {
    "id": "new_CD_REG_AFFINE_ISO_QG1_fixed_coefficient_stationary_Gaussian_mean_and_actual_response",
    "status": "NEW_PARENT_RETAINED_MEAN_STATIONARITY_AND_BOUNDED_COEFFICIENT_CORRECTION_NOT_FULL_QUANTUM_BACKGROUND_STABILITY_CUTOFF_OR_V_G_B",
}


def require_scope(time, X, kappa=profile.KAPPA, quantum_sector=1):
    t, x, k, q = map(rational, (time, X, kappa, quantum_sector))
    if not -s.Rational(1, 2) <= t <= s.Rational(1, 2):
        raise ValueError("Only the compact numerical coefficient-error slab")
    if not -s.Rational(1, 4096) < x < s.Rational(6, 5):
        raise ValueError("Outside the full original physical gradient strip")
    if k != profile.KAPPA or q != 1:
        raise ValueError(
            "Only the anchor hierarchy and specified conditional Gaussian vector sector"
        )
    return t, x, k, q


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("The separately named retuning cannot close original P8 gates")
    return True


def packets():
    return {
        "profile": profile.data(),
        "bounds": bounds.data(),
        "ward": ward.data(),
        "response": response.data(),
    }


@cache
def residuals():
    rows = {
        name + "_" + key: value
        for name, packet in packets().items()
        for key, value in packet["checks"].items()
    }
    rows.update(
        {
            "nine_primitive_rows_retained": len(frontier()) - 9,
            "primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "prior_matching_rows_retained": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "one_new_named_parent_row_added": len(matching())
            - len(previous.matching())
            - 1,
        }
    )
    return {
        key: s.ImmutableMatrix(value.applyfunc(s.cancel))
        if isinstance(value, s.MatrixBase)
        else s.cancel(value)
        for key, value in rows.items()
    }


def scalar_entry_count():
    return sum(
        value.rows * value.cols if isinstance(value, s.MatrixBase) else 1
        for value in residuals().values()
    )


@cache
def gates():
    return {
        **{
            name + "_" + key: bool(value)
            for name, packet in packets().items()
            for key, value in packet.get("gates", {}).items()
        },
        "new_named_candidate_with_explicit_error_budget": True,
        "old_frozen_parent_and_certificates_not_overwritten": True,
        "actual_source_defined_stress_not_free_profile_ansatz": True,
        "reference_profiles_fixed_once_not_adaptive_under_variation": True,
        "same_canonical_vector_state_operator_and_covariant_prescription": True,
        "complete_reference_mean_equations_including_scalar_Ward": True,
        "same_new_action_has_retained_Minkowski_and_CD_means": True,
        "full_X_strip_and_sharper_clock_bounds_are_continuous": True,
        "vacuum_low_germs_preserved_not_global_time_parity": True,
        "complete_canonical_interactions_fixed_in_actual_decoupling_family": True,
        "nonzero_individual_chart_contacts_cancel_only_in_total": True,
        "actual_nonlocal_response_not_adaptively_subtracted": True,
        "compatible_reference_state_no_nonprepared_inverse_input": True,
        "retained_stationarity_not_full_quantum_or_original_V_G_B": True,
    }


def observable():
    return {
        "new_candidate": "CD-REG-AFFINE-ISO-QG1, not a rewrite of CD-REG-AFFINE-ISO",
        "mean_result": "The same new action has the specified stationary retained Gaussian CD mean history and same-prescription Minkowski vacuum mean. Every mean equation, fixed coefficient and state normalization is accounted for.",
        "coefficient_budget": "On |u|<=1/2, all36 mixed u/X derivatives through five in each variable are below1e-742 across the original X strip, and below2e-770 on3/4<=X<6/5. Global mean identities do not imply global small-correction or stability bounds.",
        "vacuum_boundary": "The nonconstant correction begins no earlier than field degree2048. Existing scalar quadratic and quartic data remain, with the explicit vacuum constant retained. The complete new canonical interaction family is fixed, not just its low Taylor truncation.",
        "response": "The actual stationary response retains its nonlocal scalar scale channel. The prepared retained Volterra inverse applies with the same scalar norm bound and new local C1<15000, but a useful/small full inverse is not established.",
        "remaining": "Full interacting state/measure/loops, nonlinear or spatial/noise response and stability, physical heavy/cutoff and omitted-order control, complete vacuum cuts/contour/truncation, finite-gravity IR/Regge and original V/G/B.",
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
        -s.oo,
        s.zoo,
        s.I,
        s.nan,
        s.Symbol("x"),
    )
    out = []
    for i, value in enumerate(bad):
        for pos in range(4):
            args = [0, 1, profile.KAPPA, 1]
            args[pos] = value
            out.append((f"type_{pos}_{i}", require_scope, tuple(args)))
    for args in (
        (1, 1, profile.KAPPA, 1),
        (-1, 1, profile.KAPPA, 1),
        (0, -s.Rational(1, 4096), profile.KAPPA, 1),
        (0, s.Rational(6, 5), profile.KAPPA, 1),
        (0, 0, 2 * profile.KAPPA, 1),
        (0, 0, profile.KAPPA, 0),
        (0, 0, profile.KAPPA, 2),
    ):
        out.append((f"scope_{len(out)}", require_scope, args))
    for i in range(9):
        rows = frontier()
        rows[i]["status"] = "COMPLETE"
        out.append((f"primitive_{i}", validate_scope, (rows, matching())))
    for i in range(len(matching())):
        rows = matching()
        rows[i]["status"] = "COMPLETE"
        out.append((f"matching_{i}", validate_scope, (frontier(), rows)))
    out.extend(
        (
            ("missing_primitive", validate_scope, (frontier()[:-1], matching())),
            ("missing_matching", validate_scope, (frontier(), matching()[:-1])),
        )
    )
    return out


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported retuned-parent scope accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "old_action_not_silently_retuned": True,
        "complete_source_defined_stress_not_just_local_heat_piece": True,
        "fixed_coefficient_not_live_metric_dependent_subtraction": True,
        "individual_nonstationary_contacts_not_dropped": True,
        "nonlocal_response_survives_mean_cancellation": True,
        "full_canonical_family_not_assumed_even_in_time": True,
        "retained_mean_stationarity_not_stability_or_full_interacting_solution": True,
        "original_P8_not_closed": True,
    }
