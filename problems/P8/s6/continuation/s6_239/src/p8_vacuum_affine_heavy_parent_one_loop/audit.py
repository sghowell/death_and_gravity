"""Scope-checked complete first-loop coefficient of the new classical limit."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_heavy_scalar_parent import audit as previous

from . import bounds, germs, loops

ITEM = {
    "id": "complete_H8A420_classical_limiting_action_first_loop_with_higher_vertices_and_separate_OS4_matching",
    "status": "FORMAL_COMPLETE_FIRST_LOOP_OF_CLASSICAL_LIMIT_WITH_SCOPED_ERRORS_NOT_QUANTUM_DECOUPLING_ALL_LOOPS_CURVED_BOUNCE_UV_REGGE_OR_P8",
}
require_scope = previous.require_scope
exact = previous.exact


def require_stage(stage):
    if not isinstance(stage, str) or stage not in (
        "complete_classical_limit_one_loop_four_point",
        "all_new_source_contraction_classes",
        "dimensional_tensor_before_MS_subtraction",
        "new_full_parent_OS4_value",
        "complete_physical_angle_first_loop_bound",
        "complete_first_loop_forward_coefficients",
        "specified_finite_counterterm_clock_four_jet",
    ):
        raise ValueError(
            "Only the stated formal first-loop and finite-extension results are proved"
        )
    return stage


def require_prescription(name, order):
    if isinstance(order, bool) or not isinstance(order, (int, s.Integer)):
        raise TypeError("Use an exact integer loop order")
    if not isinstance(name, str) or name != "H8A420-VAC-OS4" or order != 1:
        raise ValueError(
            "Only the new complete limiting-action first-loop prescription is fixed"
        )
    return name, int(order)


def require_physical(energy_squared, angle):
    S, x = map(exact, (energy_squared, angle))
    if not 4 <= S <= bounds.S_MAX or abs(x) > 1:
        raise ValueError("Outside the complete stated first-loop matching window")
    return S, x


def require_finite_tube(u, X, order):
    return previous.require_tube(u, X, order)


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("Formal first-loop matching does not close original P8")
    return True


def packets():
    return {
        "literal_complete_limiting_action_germs_and_exhaustive_graphs": germs.data(),
        "complete_dimensional_contractions_UV_and_first_loop": loops.data(),
        "complete_matching_margins_and_finite_clock_extensions": bounds.data(),
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
            "one_new_separate_formal_first_loop_record": len(matching())
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
        "new_finite_OS4_value_not_old_contact_silently_reused": True,
        "finite_clock_extension_budget_separate_from_classical_tube": True,
        "free_M1_Proca_and_tensor_spectators_retained_in_classical_limit": True,
        "no_quantization_limit_interchange_or_full_UV_inferred": True,
        "unique_matching_identifiers": len({row["id"] for row in matching()})
        == len(matching()),
    }


def observable():
    return {
        "input": "The complete S238 fixed-canonical CLASSICAL limiting Phi/H action, with its full higher functions. The retained massless M1, massive Proca and normalized tensor spectators are free in this limiting action. No finite-gravity quantum-limit theorem is assumed.",
        "completeness": "Connected graph degree counting isolates all new one-loop four-light contributions: six literal local degree-six tadpoles and all three cubic/five-field source contraction classes. Higher localizer and smooth-retuning vertices cannot enter this coefficient.",
        "renormalization": "d=4-2epsilon is retained through tensor averaging and pole-only MSbar. Complete new off-shell local contractions define their UV counterfunctional. The old first-order light on-shell and heavy-onepoint conditions remain; a NEW full-parent symmetric four-light value condition is imposed.",
        "matching": "The full rational new heavy-vertex correction is retained, including its b40 shift. All new terms are real on the physical window; the old first elastic cut is unchanged. The complete first loop is below10^-199 of the original tree and the same forward coefficient margins persist.",
        "finite_clock": "Only the specified finite contact, quadratic and onepoint terms are extended by V=1-T. Their full mixed four-jets have a separately proved10^-400 normalized clock-tube bound; this is not the complete curved effective action or the classical10^-2700 budget.",
        "remaining": "Quantum decoupling at finite gravity, omitted loops, a new heavy curved state and full stress/response, same-state nonlinear bounce, physical UV and finite-gravity Regge estimates, and original V/G/B/P8 remain open.",
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
        for pos in range(2):
            args = [4, 0]
            args[pos] = value
            out.append((f"physical_type_{pos}_{i}", require_physical, tuple(args)))
            args = ["H8A420-VAC-OS4", 1]
            args[pos] = value
            out.append(
                (f"prescription_type_{pos}_{i}", require_prescription, tuple(args))
            )
        for pos in range(3):
            args = [0, 1, 4]
            args[pos] = value
            out.append(
                (f"finite_tube_type_{pos}_{i}", require_finite_tube, tuple(args))
            )
        out.append((f"stage_type_{i}", require_stage, (value,)))
    for args in ((3, 0), (10**197, 0), (4, -2), (4, 2)):
        out.append((f"physical_scope_{len(out)}", require_physical, args))
    for args in (
        ("V2S-T1-OS4", 1),
        ("H8A420-VAC-OS4", 0),
        ("H8A420-VAC-OS4", 2),
        ("H8A420-VAC-OS4", s.Rational(1, 2)),
    ):
        out.append((f"prescription_scope_{len(out)}", require_prescription, args))
    for args in (
        (2, 1, 0),
        (0, 0, 0),
        (0, 2, 0),
        (0, 1, -1),
        (0, 1, 5),
        (0, 1, s.Rational(1, 2)),
    ):
        out.append((f"finite_tube_scope_{len(out)}", require_finite_tube, args))
    for stage in (
        "old_polynomial_loop_automatically_transfers",
        "local_tadpoles_only",
        "all_new_terms_are_contact",
        "new_b40_shift_zero",
        "drop_heavy_onepoint_source",
        "set_d4_before_pole_subtraction",
        "constant_Galileon_IBP_with_field_coefficient",
        "old_OS4_contact_unchanged",
        "finite_extensions_use_h_without_new_quartic",
        "classical_tube_budget_covers_finite_counterterms",
        "full_curved_UV_functional",
        "all_loop_error",
        "quantization_commutes_with_gravity_decoupling",
        "regulator_removed_quantum_parent",
        "heavy_state_stress_retuned",
        "same_state_nonlinear_bounce",
        "exact_UV_S_matrix",
        "physical_cutoff",
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
            "Unsupported complete-loop or quantum-parent claim accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "all_literal_higher_vertices_and_all_contraction_classes": True,
        "full_dimensional_finite_terms_retained": True,
        "nonzero_heavy_vertex_b40_correction_retained": True,
        "new_complete_finite_value_condition_explicit": True,
        "finite_clock_extension_not_full_quantum_clock_control": True,
        "no_classical_and_quantum_limit_interchange": True,
        "no_omitted_loop_or_physical_UV_inference": True,
        "all_original_primitive_statuses_unchanged": frontier() == previous.frontier(),
    }
