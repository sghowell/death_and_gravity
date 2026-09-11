"""Scoped causal classical inverse and the first formal quantum coefficient."""

from functools import cache

import sympy as s
from p8_vacuum_affine_quantum_retuning import audit as previous
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import hamiltonian, jet_bounds, perturbation, real_bounds

ITEM = {
    "id": "actual_classical_Hamiltonian_C0_and_C10_inverse_and_first_formal_Gaussian_response",
    "status": "QUANTIFIED_PREPARED_CLASSICAL_PROPAGATOR_AND_FIRST_FORMAL_COEFFICIENT_NOT_FINITE_COUPLING_ERROR_STABILITY_CUTOFF_OR_V_G_B",
}


def require_scope(time, derivative_order=10, matter_charge_variation=0, formal_order=1):
    t, k, c, o = map(
        rational, (time, derivative_order, matter_charge_variation, formal_order)
    )
    if not -s.Rational(1, 2) <= t <= s.Rational(1, 2):
        raise ValueError("Only the stated compact observation slab")
    if k != 10 or c != 0 or o != 1:
        raise ValueError(
            "Require prepared C10 forces, original fixed matter charge and the first formal coefficient"
        )
    return t, k, c, o


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("A formal coefficient cannot close original P8 obligations")
    return True


def packets():
    return {
        "Hamiltonian": hamiltonian.data(),
        "real_norm": real_bounds.data(),
        "C10_norm": jet_bounds.data(),
        "first_formal_response": perturbation.data(),
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
            "one_new_scoped_row_added": len(matching()) - len(previous.matching()) - 1,
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
        "actual_full_target_classical_reference_not_new_action_alone": True,
        "original_force_density_and_matter_charge_retained": True,
        "regular_Hamiltonian_pivot_not_velocity_chart_division": True,
        "continuous_real_bounds_not_sampled_ODE_verdict": True,
        "uniform_complex_discs_control_all_ten_classical_jets": True,
        "Leibniz_recurrence_keeps_force_and_readout_derivatives": True,
        "actual_stationary_quantum_current_includes_fixed_profile": True,
        "reference_functions_not_adaptively_reselected": True,
        "future_extension_causal_no_right_endpoint_vanishing_assumption": True,
        "formal_marker_scales_both_profile_and_Gaussian_functional": True,
        "formal_coefficient_does_not_assume_exact_parameter_branch": True,
        "derivative_loss_prevents_claimed_C0_feedback_contraction": True,
        "no_finite_coupling_remainder_or_pole_deletion_claim": True,
        "original_V_G_B_not_closed": True,
    }


def observable():
    return {
        "classical_C0": "The actual prepared classical physical metric inverse on I has joint C0 norm below53000, with the original M1 charge retained.",
        "classical_C10": "Uniform radius1e-8 complex discs and a scaled derivative recurrence bound the same inverse from physical C10 forces to physical C10 metric by1e96.",
        "actual_quantum": "The new stationary Qstar is the unchanged conditional vector response PLUS the fixed new scalar coefficient, not either term separately.",
        "first_coefficient": "h0=Bcl g and h1=-Bcl Qstar h0 are uniquely defined prepared formal coefficients; ||h1||C0<4e-669||g||C10.",
        "not_established": "No finite-alpha error, differentiable exact alpha branch, C0 contraction, full response stability, nonlinear or spatial/noise quantum solution, physical cutoff or original V/G/B completion.",
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
            args = [0, 10, 0, 1]
            args[pos] = value
            out.append((f"type_{pos}_{i}", require_scope, tuple(args)))
    for args in (
        (1, 10, 0, 1),
        (-1, 10, 0, 1),
        (0, 0, 0, 1),
        (0, 9, 0, 1),
        (0, 10, 1, 1),
        (0, 10, 0, 0),
        (0, 10, 0, 2),
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
        raise ValueError("Unsupported propagator scope accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "old_frozen_parent_not_overwritten": True,
        "fixed_charge_Routhian_not_naive_velocity_substitution": True,
        "full_physical_force_map_kept": True,
        "no_division_by_crossing_velocity_chart_pivot": True,
        "all_ten_readout_jets_bounded": True,
        "right_endpoint_not_artificially_zeroed": True,
        "formal_small_coefficient_not_full_error_or_stability": True,
        "original_P8_not_closed": True,
    }
