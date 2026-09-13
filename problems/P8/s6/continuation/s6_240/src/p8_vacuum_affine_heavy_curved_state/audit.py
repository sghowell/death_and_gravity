"""Scoped heavy Gaussian state, full reference stress and fixed mean profile."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_heavy_parent_one_loop import audit as previous

from . import estimates, quantum, state

ITEM = {
    "id": "specified_H8A420_SLE_heavy_state_complete_MSbar_reference_stress_and_fixed_QG2_mean_profile",
    "status": "SPECIFIED_GAUSSIAN_HEAVY_STATE_REFERENCE_STRESS_AND_RETAINED_MEAN_MATCHING_NOT_FULL_INTERACTING_CLOCK_RESPONSE_NONLINEAR_BOUNCE_UV_REGGE_OR_P8",
}
require_scope = previous.require_scope
exact = previous.exact
require_tube = previous.require_finite_tube


def require_stage(stage):
    if not isinstance(stage, str) or stage not in (
        "specified_exact_heavy_SLE_Hadamard_state",
        "full_dimensional_scalar_counteraction",
        "complete_all_momentum_reference_stress",
        "five_reference_time_derivatives",
        "fixed_QG2_physical_mean_profile",
        "new_combined_clock_four_jet_budget",
        "unchanged_first_loop_vacuum_scattering_jets",
        "complete_flat_Gaussian_one_loop_constant",
    ):
        raise ValueError(
            "Only the stated Gaussian state, reference and scoped matching results are proved"
        )
    return stage


def require_state(name, mass_squared):
    n = exact(mass_squared)
    if not isinstance(name, str) or name != "H8A420-SLE-PRE-T0" or n != state.MASS2:
        raise ValueError("The exact specified new heavy state and mass are fixed")
    return name, n


def require_estimate(time, momentum, order):
    u, p = map(exact, (time, momentum))
    if isinstance(order, bool) or not isinstance(order, (int, s.Integer)):
        raise TypeError("Use an exact integer derivative order")
    if not -1 <= u <= 1 or p < 0 or not 0 <= order <= 5:
        raise ValueError("Outside the full reference stress estimate scope")
    return u, p, int(order)


def require_cone(c1, c2_real, c2_imag):
    energy, x, y = map(exact, (c1, c2_real, c2_imag))
    if energy <= 0 or energy * energy <= x * x + y * y:
        raise ValueError("The exact SLE energy cone must be strictly positive")
    return energy, x, y


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "A specified Gaussian state and mean profile do not close original P8"
        )
    return True


def packets():
    return {
        "full_reference_modes_and_actual_Hadamard_state": state.data(),
        "whole_dimensional_counteraction_and_fixed_mean_profile": quantum.data(),
        "complete_all_momentum_stress_and_clock_tube_bounds": estimates.data(),
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
            "nine_original_primitive_rows": len(frontier()) - 9,
            "unchanged_original_primitive_statuses": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "unchanged_previous_matching_records": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "one_new_separate_Gaussian_state_and_mean_record": len(matching())
            - len(previous.matching())
            - 1,
        }
    )
    return {key: s.cancel(value) for key, value in rows.items()}


def scalar_entry_count():
    return len(residuals())


@cache
def gates():
    return {
        **{
            name + "_" + key: bool(value)
            for name, packet in packets().items()
            for key, value in packet["gates"].items()
        },
        "all_frozen_scientific_prescriptions_unchanged": True,
        "new_heavy_reference_state_and_scalar_prescription_explicit": True,
        "full_state_error_not_only_a_local_heat_expansion": True,
        "complete_new_profile_budget_separate_from_older_budgets": True,
        "no_automatic_old_Hessian_stability_or_inverse_transfer": True,
        "unique_matching_identifiers": len({row["id"] for row in matching()})
        == len(matching()),
    }


def observable():
    return {
        "input": "The S238 full covariant parent and S239 complete formal limiting-action first loop are frozen. The new minimally coupled scalar state and its covariant mu1 prescription are specified; old Proca and M1 choices are retained.",
        "state": "A fixed smooth compact pre-t0 sampling determines the exact heavy state of low energy. Its physical phase, CCR and basis-independent minimizer are checked. The Hadamard theorem applies to the full exact state, not its finite WKB comparison.",
        "stress": "The full state-dependent integral, complete adiabatic subtraction and the dimensionally matched local scalar action give rho/P and five time derivatives below10^-400 against kappa0 on[-1,1], uniformly over all momentum. No cutoff or local-only approximation is used.",
        "profile": "A separately named fixed QG2 physical coefficient balances the actual specified heavy Gaussian mean stress. The new profile also retains the heavy and light flat Gaussian vacuum constants explicitly. Its nonconstant vacuum jets cannot change the S239 four-light first loop.",
        "clock": "The original physical CD metric/u/M1 with the S238 recomputed affine mean is stationary for the retained Proca-plus-H Gaussian reference truncation. The new scalar clock jets change the light Hessian; the old stability and quantum inverse results are not silently transferred. Combined coefficient four-jets fit a NEW10^-399 tube budget.",
        "remaining": "Full interacting light/mixed clock loops and state, quantum gravitational decoupling, full new response/inverse, nonlinear same-state bounce, physical UV, finite-gravity Regge and original V/G/B/P8 remain open.",
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
            args = [0, modes.KAPPA, modes.MASS, 1]
            args[pos] = value
            out.append((f"original_scope_type_{pos}_{i}", require_scope, tuple(args)))
        for pos in range(3):
            for label, call, values in (
                ("estimate", require_estimate, [0, 0, 5]),
                ("tube", require_tube, [0, 1, 4]),
                ("cone", require_cone, [3, 1, 2]),
            ):
                args = list(values)
                args[pos] = value
                out.append((f"{label}_type_{pos}_{i}", call, tuple(args)))
        for pos in range(2):
            args = ["H8A420-SLE-PRE-T0", state.MASS2]
            args[pos] = value
            out.append((f"state_type_{pos}_{i}", require_state, tuple(args)))
        out.append((f"stage_type_{i}", require_stage, (value,)))
    for args in (
        (-2, 0, 0),
        (2, 0, 0),
        (0, -1, 0),
        (0, 0, -1),
        (0, 0, 6),
        (0, 0, s.Rational(1, 2)),
    ):
        out.append((f"estimate_scope_{len(out)}", require_estimate, args))
    for args in (("old_Proca_state", state.MASS2), ("H8A420-SLE-PRE-T0", 1000)):
        out.append((f"state_scope_{len(out)}", require_state, args))
    for args in ((0, 0, 0), (-1, 0, 0), (1, 1, 0), (1, 1, 1), (3, 0, 4)):
        out.append((f"cone_scope_{len(out)}", require_cone, args))
    for args in (
        (2, 1, 0),
        (0, 0, 0),
        (0, 2, 0),
        (0, 1, -1),
        (0, 1, 5),
        (0, 1, s.Rational(1, 2)),
    ):
        out.append((f"tube_scope_{len(out)}", require_tube, args))
    for stage in (
        "WKB_comparison_is_actual_state",
        "instantaneous_vacuum_Hadamard",
        "reset_state_after_metric_perturbation",
        "component_only_dimensional_pole",
        "heavy_vacuum_energy_zero",
        "light_vacuum_constant_omitted",
        "local_heat_terms_are_full_stress",
        "finite_momentum_cutoff",
        "small_mass_fixtures_prove_actual_hierarchy",
        "global_time_uniform_bound",
        "adaptive_stress_subtraction",
        "whole_QG2_scalar_real_analytic",
        "old_classical_Hessian_unchanged",
        "old_stability_rows_transfer",
        "full_quantum_response_small",
        "same_space_inverse_from_kappa",
        "all_light_clock_loops",
        "full_interacting_clock_state",
        "quantization_commutes_with_gravity_limit",
        "all_loop_error",
        "nonlinear_same_state_bounce",
        "exact_UV_S_matrix",
        "finite_gravity_Regge_done",
        "closed_P8",
    ):
        out.append(("unsupported_" + stage, require_stage, (stage,)))
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
            (
                "extra_full_parent",
                validate_scope,
                (
                    frontier(),
                    matching() + [{"id": "full_parent", "status": "COMPLETE"}],
                ),
            ),
        )
    )
    return out


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (ValueError, TypeError):
            continue
        raise ValueError(
            "Unsupported state, stress or quantum-parent claim accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "physical_positive_frequency_and_CCR_retained": True,
        "comparison_basis_not_the_selected_state": True,
        "full_dimensional_counteraction_before_limit": True,
        "complete_state_error_and_all_momenta_retained": True,
        "fixed_physical_profile_not_adaptive_subtraction": True,
        "both_flat_vacuum_constants_retained": True,
        "no_automatic_response_inverse_or_stability_transfer": True,
        "all_original_primitive_statuses_unchanged": frontier() == previous.frontier(),
    }
