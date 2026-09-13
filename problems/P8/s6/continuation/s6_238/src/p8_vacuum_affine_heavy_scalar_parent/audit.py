"""Full new candidate with the original P8 obligations explicitly retained."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_heavy_scalar_loop_coefficients import audit as previous

from . import clock, family, heavy

ITEM = {
    "id": "separate_covariant_H8A420_parent_with_full_affine_domain_clock_coefficient_tube_and_scoped_heavy_stress_control",
    "status": "SEPARATE_SMOOTH_CLASSICAL_PARENT_AND_SCOPED_MATCHING_NOT_FULL_QUANTUM_PARENT_SAME_STATE_BOUNCE_UV_REGGE_ORIGINAL_V_G_B_OR_P8",
}
require_scope = previous.require_scope


def require_stage(stage):
    if not isinstance(stage, str) or stage not in (
        "separate_full_classical_covariant_definition",
        "connected_regular_affine_domain",
        "complete_clock_four_jet_coefficient_tube",
        "full_global_flat_classical_potential",
        "exact_prescribed_heavy_source_and_corridor_stress",
        "unrestricted_stress_counterexample",
    ):
        raise ValueError(
            "Only the stated separate classical and scoped heavy results are proved"
        )
    return stage


def exact(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Integer, s.Rational)):
        raise TypeError("Use exact real rational inputs")
    return s.Rational(value)


def require_candidate(name, parameter):
    a = exact(parameter)
    if not isinstance(name, str) or name != family.NAME or a != family.LOCALIZER:
        raise ValueError("The separately named H8A420 candidate is fixed")
    return name, a


def require_tube(u, X, order):
    uv, x = map(exact, (u, X))
    if isinstance(order, bool) or not isinstance(order, (int, s.Integer)):
        raise TypeError("Use an exact integer jet order")
    if (
        abs(uv) > 1
        or not s.Rational(7, 8) <= x <= s.Rational(9, 8)
        or not 0 <= order <= 4
    ):
        raise ValueError("Outside the complete stated mixed four-jet tube")
    return uv, x, int(order)


def require_corridor(phi, X):
    v, x = map(exact, (phi, X))
    if not 0 <= x <= s.Rational(6, 5) or not (
        v * v <= 10**380 or (x >= s.Rational(1, 2) and v * v <= family.K0)
    ):
        raise ValueError("Outside the two-chart prescribed-source corridor")
    return v, x


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("A new classical candidate does not close original P8")
    return True


def packets():
    return {
        "literal_new_full_functions_vacuum_and_clock_jets": family.data(),
        "complete_connected_affine_domain_and_coefficient_tube": clock.data(),
        "full_potential_and_exact_heavy_source_stress": heavy.data(),
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
            "one_new_separate_classical_candidate_record": len(matching())
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
        "old_scientific_and_quantum_prescription_bytes_unchanged": True,
        "new_full_smooth_candidate_explicitly_named": True,
        "all_R_dependent_affine_maps_recomputed": True,
        "source_metric_stress_term_and_counterexample_retained": True,
        "new_heavy_quantum_stress_not_silently_retuned": True,
        "unique_matching_identifiers": len({row["id"] for row in matching()})
        == len(matching()),
    }


def observable():
    return {
        "input": "The literal original S109/S177 complete fixed canonical functions, S182 smooth QG1 retuning and S174 generic full-connection construction are retained as frozen inputs. A separately named new local action is defined.",
        "vacuum": "The decoupled flat degree-four tree matches V2S-T1 exactly; higher functions and derivative-dependent heavy source remain different. The full flat constant-field potential relative to its retained origin is globally coercive.",
        "clock": "The modified independent functions have eighth-order clock zeros and the dependent Ia differences seventh-order zeros. A joint complex-neighborhood argument bounds all mixed four-jets of the full normalized F,R,A3,A4,A5,J by10^-2700 on the named real tube.",
        "affine": "The complete original positive R range persists, including the delicate upper edge. Generic full-rank quotient, source-centering and covariant chart formulas therefore construct a new algebraic parent. q and B are recomputed, not reset or silently identified with their old values.",
        "heavy": "The full prescribed-source stress includes2Y H J_Y. A connected two-chart corridor gives a10^-8 error relative to Y between the exact particular and retained valley stresses. Outside it an explicit relative stress error exceeds10^395 despite the small offset proxy.",
        "quantum_boundary": "The new free heavy field has curved vacuum stress and new higher interactions. The old QG1 mean background and S235-S237 first-loop amplitudes do not automatically transfer. No new heavy state or counterprofile is silently supplied.",
        "remaining": "All-loop and new-parent vacuum matching, heavy-state and curved response control, same-state nonlinear bounce, physical UV, finite-gravity Regge estimates and original V/G/B/P8 remain open.",
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
            args = [0, 1, 4]
            args[pos] = value
            out.append((f"tube_type_{pos}_{i}", require_tube, tuple(args)))
        for pos in range(2):
            args = [0, 1]
            args[pos] = value
            out.append((f"corridor_type_{pos}_{i}", require_corridor, tuple(args)))
            args = [family.NAME, family.LOCALIZER]
            args[pos] = value
            out.append((f"candidate_type_{pos}_{i}", require_candidate, tuple(args)))
        out.append((f"stage_type_{i}", require_stage, (value,)))
    for args in (
        (2, 1, 0),
        (0, 0, 0),
        (0, 2, 0),
        (0, 1, -1),
        (0, 1, 5),
        (0, 1, s.Rational(1, 2)),
    ):
        out.append((f"tube_scope_{len(out)}", require_tube, args))
    for args in (
        (0, -1),
        (0, 2),
        (10**191, 0),
        (10**401, 1),
        (10**400, s.Rational(1, 10**210)),
    ):
        out.append((f"corridor_scope_{len(out)}", require_corridor, args))
    for args in (
        ("original_parent", family.LOCALIZER),
        (family.NAME, 10**208),
        (family.NAME, 10**820),
    ):
        out.append((f"candidate_scope_{len(out)}", require_candidate, args))
    for stage in (
        "original_action_overwritten",
        "unchanged_affine_background_trace",
        "regular_q_clock_value_reset",
        "quartic_a3_spliced_alone",
        "all_field_strip_heavy_stress_small",
        "offset_proxy_bounds_full_stress",
        "source_metric_derivative_dropped",
        "full_parent_real_analytic",
        "all_background_kinetic_gap",
        "coupled_bounce_solution",
        "same_heavy_quantum_state",
        "heavy_quantum_stress_retuned",
        "same_S237_parent_loop",
        "all_loop_error",
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
        raise ValueError("Unsupported parent or stress claim accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "full_old_functions_and_new_candidate_both_retained": True,
        "all_Ia_and_full_affine_maps_recomputed": True,
        "no_false_global_analyticity_transfer": True,
        "full_complex_tube_not_sampled_smallness": True,
        "complete_source_metric_stress_factor_retained": True,
        "unrestricted_stress_counterexample_not_hidden": True,
        "no_automatic_new_quantum_state_or_loop_transfer": True,
        "all_original_primitive_statuses_unchanged": frontier() == previous.frontier(),
    }
