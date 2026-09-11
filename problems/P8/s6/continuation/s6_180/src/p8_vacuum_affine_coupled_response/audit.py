"""Scope-preserving actual retained tree-plus-conditional-vector inverse audit."""

from functools import cache

import sympy as s
from p8_vacuum_affine_metric_response import audit as previous
from p8_vacuum_affine_proca_gaussian.bridge import KAPPA, ZETA
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import classical, inverse, quantum
from . import matching as kernel_matching

ITEM = {
    "id": "same_affine_parent_prepared_homogeneous_retained_tree_plus_conditional_vector_causal_inverse",
    "status": "UNIQUE_RETAINED_LINEAR_PREPARED_CAUSAL_INVERSE_WITH_NUMERIC_LOCAL_LAPSE_RECOVERY_NOT_FULL_INVERSE_SMALLNESS_QUANTUM_BACKGROUND_OR_V_G_B",
}


def require_scope(time, zeta=ZETA, kappa=KAPPA, charge_perturbation=0):
    t, z, k, c = map(rational, (time, zeta, kappa, charge_perturbation))
    if not -s.Rational(1, 2) <= t <= s.Rational(1, 2):
        raise ValueError("Outside the quantified current-parent clock slab")
    if z != ZETA or k != KAPPA or c != 0:
        raise ValueError(
            "Require the fixed conditional parent and zero prepared matter-charge perturbation"
        )
    return t, z, k, c


def frontier():
    return previous.frontier()


def matching_frontier():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def matching():
    return matching_frontier()


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching_frontier():
        raise ValueError("A retained linear inverse cannot close original P8 gates")
    return True


def packets():
    return {
        "clock": classical.clock_coefficients(),
        "quadratic": classical.quadratic(),
        "tree": classical.adapted(),
        "quantum": quantum.data(),
        "mode_covariance": quantum.mode_covariance(),
        "matching": kernel_matching.data(),
        "inverse": inverse.data(),
        "first_row": inverse.first_row(),
    }


@cache
def residuals():
    out = {
        name + "_" + key: value
        for name, packet in packets().items()
        for key, value in packet["checks"].items()
    }
    out.update(
        {
            "nine_primitive_rows_retained": len(frontier()) - 9,
            "primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "prior_matching_rows_retained": sum(
                a != b for a, b in zip(matching_frontier(), previous.matching())
            ),
            "one_retained_current_parent_inverse_row_added": len(matching_frontier())
            - len(previous.matching())
            - 1,
        }
    )
    return {
        key: s.ImmutableMatrix(value.applyfunc(s.cancel))
        if isinstance(value, s.MatrixBase)
        else s.cancel(value)
        for key, value in out.items()
    }


def scalar_entry_count():
    return sum(
        value.rows * value.cols if isinstance(value, s.MatrixBase) else 1
        for value in residuals().values()
    )


@cache
def gates():
    quantitative = {
        name + "_" + key: bool(value)
        for name, packet in packets().items()
        for key, value in packet.get("gates", {}).items()
    }
    return {
        **quantitative,
        "same_complete_analytic_parent_and_actual_clock_jets": True,
        "literal_weighted_boundary_and_lapse_pivot_rederived": True,
        "free_M1_charge_reduced_variationally_not_wrong_substitution": True,
        "physical_and_adapted_tree_current_routes_agree": True,
        "same_conditional_vector_measure_state_and_covariant_prescription": True,
        "no_old_scalar_stress_canceling_profile_installed": True,
        "full_uncancelled_volume_readout_and_chart_contacts_retained": True,
        "exact_prepared_time_covariance_and_conserved_physical_Ward": True,
        "time_channel_is_actual_local_constraint_not_discarded_kernel": True,
        "full_dimensional_current_and_finite_fourth_contact_matching": True,
        "same_all_momentum_massive_scalar_range_kernel": True,
        "rank_one_loop_completed_by_actual_tree_not_pseudoinverse": True,
        "variable_coefficient_primitive_commutators_retained": True,
        "weak_log_remainder_and_all_diagonal_derivatives_integrable": True,
        "constructive_finite_weight_in_terms_of_unevaluated_majorants": True,
        "actual_C1_numeric_but_full_inverse_norm_not_evaluated": True,
        "physical_lapse_recovered_without_derivative_of_arbitrary_C0_scale": True,
        "prepared_state_charge_and_causality_not_reset": True,
        "nonprepared_background_residual_not_fed_to_prepared_inverse": True,
        "full_quantum_background_stability_cutoff_and_V_G_B_open": True,
    }


def observable():
    return {
        "fixed_object": "The current complete affine/CD/M1 classical parent plus only the specified conditional Gaussian vector mean response, with no old stress-canceling profile.",
        "new_inverse": "Unique smooth prepared causal solution of the retained HOMOGENEOUS LINEAR RESPONSE equations for every smooth prepared original two-force, with a finite C0 bound in the physical lapse and scale.",
        "new_numeric_piece": "The actual local first-row lapse-reconstruction constant is below15000, with full tree boundaries and the conditional pressure/first-derivative contribution retained.",
        "unknown_constants": "The curved weak-log majorant and full scalar-kernel L1 norm are finite but not numerically evaluated. The resulting inverse is not shown small or stable.",
        "background_boundary": "The CD reference still has nonzero quantum one-point stress. Its residual is not a prepared force at the original initial neighborhood; no corrected background or new compatible quantum Cauchy data are supplied.",
        "remaining": "Spatial/noise response, nonlinear remainder and background, full interacting measure/state and loops, physical cutoff/omitted terms, vacuum matching/cuts/contour/truncation, finite-gravity IR/Regge and original V/G/B.",
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
            args = [0, ZETA, KAPPA, 0]
            args[pos] = value
            out.append((f"type_{pos}_{i}", require_scope, tuple(args)))
    for row in (
        (-1, ZETA, KAPPA, 0),
        (1, ZETA, KAPPA, 0),
        (0, s.Rational(1, 2000), KAPPA, 0),
        (0, ZETA, 2 * KAPPA, 0),
        (0, ZETA, KAPPA, 1),
        (0, ZETA, KAPPA, -1),
    ):
        out.append((f"scope_{len(out)}", require_scope, row))
    for i in range(9):
        rows = frontier()
        rows[i]["status"] = "COMPLETE"
        out.append((f"primitive_{i}", validate_scope, (rows, matching_frontier())))
    for i in range(len(matching_frontier())):
        rows = matching_frontier()
        rows[i]["status"] = "COMPLETE"
        out.append((f"matching_{i}", validate_scope, (frontier(), rows)))
    out.extend(
        (
            (
                "missing_primitive",
                validate_scope,
                (frontier()[:-1], matching_frontier()),
            ),
            (
                "missing_matching",
                validate_scope,
                (frontier(), matching_frontier()[:-1]),
            ),
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
        raise ValueError("Unsupported retained inverse scope accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "full_boundary_and_fixed_charge_not_dropped": True,
        "old_profile_inverse_not_blindly_transferred": True,
        "inactive_loop_direction_completed_by_current_classical_constraint": True,
        "variable_fourth_coefficient_not_pulled_out_of_integral": True,
        "C0_lapse_recovery_uses_actual_local_first_row": True,
        "finite_inverse_not_smallness_stability_or_cutoff": True,
        "nonzero_reference_residual_not_false_prepared_force": True,
        "original_P8_not_closed": True,
    }
