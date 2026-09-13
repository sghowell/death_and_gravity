"""Separate first-loop observables, without upgrading any original P8 primitive."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_heavy_scalar_tree_matching import audit as previous

from . import gaussian, self_energy, width

ITEM = {
    "id": "separate_V2S_T1_exact_heavy_integration_complete_first_light_self_energy_complex_disk_bound_and_leading_heavy_decay_coefficient",
    "status": "SEPARATE_FINITE_REGULATOR_IDENTITY_AND_FIRST_LOOP_INPUT_NOT_ALL_LOOP_FOUR_POINT_ORIGINAL_PARENT_BOUNCE_UV_OR_P8",
}
require_scope = previous.require_scope


def require_stage(stage):
    if not isinstance(stage, str) or stage not in (
        "exact_finite_regulator_heavy_integration",
        "complete_nonlocal_one_light_loop_generator",
        "first_light_on_shell_self_energy",
        "one_loop_truncated_complex_disk_comparison",
        "first_heavy_absorptive_decay_coefficient",
    ):
        raise ValueError(
            "Only the specified separate-model finite-regulator and first-loop input is proved"
        )
    return stage


def require_observable(name, invariant, order):
    if not isinstance(name, str) or name != gaussian.NAME:
        raise ValueError("Only the separate V2S-T1 model is evaluated")
    if isinstance(order, bool) or not isinstance(order, (int, s.Integer)) or order != 1:
        raise ValueError("Only the first loop coefficient is evaluated")
    if isinstance(invariant, bool) or not isinstance(invariant, (int, s.Expr)):
        raise TypeError(
            "The invariant must be exact with rational real and imaginary parts"
        )
    value = s.sympify(invariant)
    real, imag = value.as_real_imag()
    if not isinstance(real, s.Rational) or not isinstance(imag, s.Rational):
        raise TypeError(
            "The invariant must be exact with rational real and imaginary parts"
        )
    if (real - 1) ** 2 + imag**2 > self_energy.RADIUS**2:
        raise ValueError("Outside the proved complex invariant disk")
    return name, value, 1


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "First-loop input does not close the full amplitude or original P8"
        )
    return True


def packets():
    return {
        "exact_heavy_integration_and_complete_light_loop_generator": gaussian.data(),
        "complete_first_light_self_energy_and_complex_disk_bound": self_energy.data(),
        "first_heavy_absorptive_coefficient_and_decay_normalization": width.data(),
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
            "unchanged_primitive_rows": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "unchanged_previous_matching_rows": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "one_new_first_loop_input": len(matching()) - len(previous.matching()) - 1,
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
        "original_frozen_model_and_target_unchanged": True,
        "one_point_condition_and_full_mixed_loop_retained": True,
        "first_loop_subtraction_prescription_explicit": True,
        "complex_disk_proof_not_finite_sample": True,
        "one_loop_pole_and_width_not_exact_quantum_claim": True,
        "unique_matching_identifiers": len({row["id"] for row in matching()})
        == len(matching()),
    }


def observable():
    return {
        "model": "V2S-T1 only, with its frozen canonical two-scalar tree action and parameters unchanged.",
        "integration": "Exact finite-regulator Euclidean heavy integration gives a bounded nonlocal light action. The complete light Hessian retains mixed loops despite the light-independent flat heavy determinant. Its formal first-loop generator is not evaluated four-point matching.",
        "first_light_loop": "With an explicit one-point-zero and light on-shell subtraction prescription, the complete first light self-energy includes the quartic tadpole and heavy-light bubble. The latter has a uniform complex logarithm remainder, giving a one-loop-truncated propagator ratio error below10^-209 on abs(s-1)<=10^196.",
        "renormalization_boundary": "The first coefficients satisfy abs(Pi_MS(1))<10^-7 and0<Pi_MS'(1)<10^-207. The subtracted pole location and unit residue are renormalization conditions at this order, not a full quantum LSZ theorem.",
        "heavy_coefficient": "The correctly normalized first heavy absorptive coefficient is positive and gives10^-208<Gamma_first/M_H<10^-207, with the formal outgoing second-sheet sign. No exact width or stable full quantum atom is asserted.",
        "remaining": "All-loop control, evaluated complete real four-point matching, higher light cuts, the continuum quantum construction, original common-parent bounce, finite-gravity Regge control and original V/G/B/P8 remain open.",
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
            out.append((f"target_scope_type_{pos}_{i}", require_scope, tuple(args)))
        out.append((f"stage_type_{i}", require_stage, (value,)))
        out.append((f"model_type_{i}", require_observable, (value, 1, 1)))
        out.append(
            (f"loop_order_type_{i}", require_observable, (gaussian.NAME, 1, value))
        )
        if value != s.I:
            out.append(
                (f"invariant_type_{i}", require_observable, (gaussian.NAME, value, 1))
            )
    for args in (
        (1, modes.KAPPA, modes.MASS, 1),
        (-1, modes.KAPPA, modes.MASS, 1),
        (0, 2 * modes.KAPPA, modes.MASS, 1),
        (0, modes.KAPPA, 999, 1),
        (0, modes.KAPPA, modes.MASS, 2),
    ):
        out.append((f"target_scope_{len(out)}", require_scope, args))
    for args in (
        ("original_affine_parent", 1, 1),
        ("V2S-T2", 1, 1),
        ("", 1, 1),
        (gaussian.NAME, 1, 0),
        (gaussian.NAME, 1, 2),
        (gaussian.NAME, 1, -1),
        (gaussian.NAME, 1 + 2 * self_energy.RADIUS, 1),
        (gaussian.NAME, 1 + 2 * s.I * self_energy.RADIUS, 1),
        (gaussian.NAME, s.sqrt(2), 1),
    ):
        out.append((f"observable_domain_{len(out)}", require_observable, args))
    for stage in (
        "all_loop_error_below_10_minus_209",
        "exact_LSZ",
        "full_quantum_UV",
        "original_parent_replaced",
        "full_real_four_point_matching",
        "physical_cutoff",
        "common_parent_bounce",
        "curved_heavy_determinant_zero",
        "mixed_loops_absent",
        "heavy_tadpole_discarded",
        "stable_exact_heavy_atom",
        "exact_resonance_width",
        "first_sheet_resonance_pole",
        "retarded_one_copy_action",
        "no_lower_higher_loop_light_cuts",
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
                "extra_quantum_parent",
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
        raise ValueError("Unsupported first-loop closure accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "finite_regulator_identity_not_continuum_quantum_theory": True,
        "mixed_loops_not_removed_with_heavy_determinant": True,
        "one_point_counterterm_cancellation_explicit": True,
        "one_loop_on_shell_condition_not_exact_LSZ": True,
        "first_width_not_exact_stable_atom_or_resummation": True,
        "two_point_error_not_full_four_point_error": True,
        "original_bounce_and_Regge_not_supplied": True,
        "all_original_primitive_statuses_unchanged": frontier() == previous.frontier(),
    }
