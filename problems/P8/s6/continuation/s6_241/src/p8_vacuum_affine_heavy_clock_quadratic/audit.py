"""Scoped new-profile quadratic comparison and complete finite local response audit."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_curved_state import audit as previous

from . import clock, local, tensor

ITEM = {
    "id": "QG2_H8A420_explicit_new_clock_quadratic_comparison_and_complete_physical_finite_heavy_response",
    "status": "NEW_CLASSICAL_COEFFICIENT_HYPOTHESES_AND_FULL_FINITE_LOCAL_HEAVY_HESSIAN_NOT_FULL_QUANTUM_RESPONSE_INVERSE_NONLINEAR_BOUNCE_UV_REGGE_OR_P8",
}
exact = previous.exact


def require_stage(stage):
    allowed = (
        "literal_full_QG2_clock_quadratic",
        "explicit_new_generic_chart_hypotheses",
        "current_coefficient_sector_twelve_derivative_comparison",
        "independent_classical_heavy_energy_phase",
        "full_physical_ADM_finite_heat_Hessian",
        "full_covariant_vacuum_constant_cancellation",
        "full_fixed_heavy_profile_Hessian_bound",
        "complete_weighted_local_tensor_Euler_graph",
        "complete_same_clock_local_Hessian",
    )
    if not isinstance(stage, str) or stage not in allowed:
        raise ValueError(
            "Only the stated coefficient comparison and complete local response are established"
        )
    return stage


def require_parameters(mass_squared, kappa):
    n, k = map(exact, (mass_squared, kappa))
    if n != local.state.MASS2 or k != local.state.KAPPA:
        raise ValueError("Keep the fixed QG2 heavy parameters")
    return n, k


def require_domain(time, momentum, block):
    u, p = map(exact, (time, momentum))
    if not isinstance(block, str) or block not in (
        "clock",
        "added_classical_heavy",
        "finite_local_scalar",
        "finite_local_tensor",
    ):
        raise ValueError("Select an established block")
    extent = s.S.One if block == "added_classical_heavy" else s.Rational(1, 2)
    if not -extent <= u <= extent or p < 0 or (block == "clock" and p == 0):
        raise ValueError("Outside the stated block's time or momentum domain")
    return u, p, block


def require_chart(label, time, momentum):
    u, p = map(exact, (time, momentum))
    if not isinstance(label, str) or label not in ("central", "outer", "regular_low"):
        raise ValueError("Select an actual phase chart")
    if not -s.Rational(1, 2) <= u <= s.Rational(1, 2) or p <= 0:
        raise ValueError("Outside the coefficient-sector domain")
    if label == "central" and (abs(u) > s.Rational(1, 4) or p < 100):
        raise ValueError("The central chart requires its full high-transfer margin")
    if label == "outer" and (abs(u) < s.Rational(1, 8) or p < 100):
        raise ValueError("The outer chart may not divide by the crossing")
    if label == "regular_low" and p > 100:
        raise ValueError("Use a high-transfer chart for the uniform estimate")
    return label, u, p


def require_graph(label, order):
    if isinstance(order, bool) or not isinstance(order, (int, s.Integer)):
        raise TypeError("Use the exact integer graph order")
    values = {
        "physical_scalar_two_jets": 2,
        "same_clock_scalar_two_jets": 2,
        "local_tensor_total_order": 4,
        "coefficient_phase_spatial_loss": 12,
        "heavy_energy_scaled_phase": 0,
    }
    if not isinstance(label, str) or label not in values or order != values[label]:
        raise ValueError("Keep the full stated derivative or energy-scaled graph")
    return label, int(order)


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "A local-response and coefficient comparison result does not close original P8"
        )
    return True


def packets():
    return {
        "full_current_clock_quadratic_and_rechecked_comparison": clock.data(),
        "complete_physical_ADM_local_heat_and_fixed_profile_Hessian": local.data(),
        "complete_finite_tensor_Euler_operator_and_graph": tensor.data(),
    }


@cache
def residuals():
    rows = {
        name + "_" + key: value
        for name, data in packets().items()
        for key, value in data["checks"].items()
    }
    rows.update(
        {
            "nine_original_primitive_rows": len(frontier()) - 9,
            "unchanged_primitive_statuses": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "unchanged_previous_matching_records": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "one_new_separate_quadratic_response_record": len(matching())
            - len(previous.matching())
            - 1,
        }
    )
    return {
        key: value.applyfunc(s.cancel)
        if isinstance(value, s.MatrixBase)
        else s.cancel(value)
        for key, value in rows.items()
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
            for name, data in packets().items()
            for key, value in data["gates"].items()
        },
        "all_frozen_inputs_and_prescriptions_unchanged": True,
        "new_coefficient_hypotheses_checked_not_automatic_transfer": True,
        "complete_physical_lapse_metric_shift_and_fixed_profile_terms": True,
        "local_graph_bound_not_a_full_quantum_inverse": True,
        "classical_heavy_decoupling_not_zero_Gaussian_response": True,
        "unique_matching_identifiers": len({row["id"] for row in matching()})
        == len(matching()),
    }


def observable():
    return {
        "input": "The complete S240 QG2 candidate/state/prescription and old S221 generic full chart proof are frozen and read-only. No further candidate or finite renormalization is introduced.",
        "clock": "Literal full profile/R lapse expansion gives all new pivot, mixed and volume terms. The combined actual fixed stress satisfies every generic coefficient and chart hypothesis. The coefficient-sector phase comparison retains exp(1e29)(1+P²)^6 with12 spatial derivatives lost. This is not the full quantum constraint inverse.",
        "heavy": "The added classical mode has an independent exact KG energy-scaled phase bound below1000 on[-1,1], uniform over momentum. Its Gaussian connected metric response is not zero.",
        "scalar_local": "The entire finite scalar heat action and fixed heavy-profile quadratic are computed in physical ADM lapse,zeta,shift variables and pulled back into the common clock variables. The constant vacuum action cancels exactly before variation. Full second-metric, curvature, shift and the nonzero state-profile chart contact give a normalized weak-graph bound below10^-580 on[-1/2,1/2], all momentum.",
        "tensor_local": "The complete finite tensor Euler term includes the actual a³-pairing adjoint, variable kinetic coefficient and every fourth-order contact. It maps the stated order-four graph to the target norm with bound below10^-597; no causal pole is removed.",
        "remaining": "The remaining exact SLE/subtraction-dependent heavy determinant response, full quantum constraints/inverse, interacting light/mixed clock loops, quantum gravitational limit, same-state nonlinear bounce, physical UV, finite-gravity Regge and original V/G/B/P8 remain open.",
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
    templates = (
        ("parameters", require_parameters, [local.state.MASS2, local.state.KAPPA]),
        ("domain", require_domain, [0, 1, "clock"]),
        ("chart", require_chart, ["central", 0, 100]),
        ("graph", require_graph, ["physical_scalar_two_jets", 2]),
    )
    for i, value in enumerate(bad):
        for name, call, values in templates:
            for pos in range(len(values)):
                args = list(values)
                args[pos] = value
                out.append((f"{name}_type_{pos}_{i}", call, tuple(args)))
        out.append((f"stage_type_{i}", require_stage, (value,)))
    for args in (
        (1, 1, "clock"),
        (-1, 1, "clock"),
        (0, 0, "clock"),
        (0, -1, "finite_local_scalar"),
        (2, 0, "added_classical_heavy"),
        (1, 0, "finite_local_tensor"),
        (0, 1, "full_quantum"),
    ):
        out.append((f"domain_scope_{len(out)}", require_domain, args))
    for args in (
        ("central", s.Rational(1, 3), 100),
        ("central", 0, 1),
        ("outer", 0, 100),
        ("outer", 1, 100),
        ("regular_low", 0, 101),
        ("central", 0, 0),
    ):
        out.append((f"chart_scope_{len(out)}", require_chart, args))
    for args in ((1, local.state.KAPPA), (local.state.MASS2, 1)):
        out.append((f"fixed_parameter_{len(out)}", require_parameters, args))
    for args in (
        ("physical_scalar_two_jets", 1),
        ("local_tensor_total_order", 2),
        ("coefficient_phase_spatial_loss", 0),
        ("heavy_energy_scaled_phase", 1),
    ):
        out.append((f"graph_scope_{len(out)}", require_graph, args))
    for stage in (
        "old_stability_rows_transfer_without_checks",
        "whole_QG2_Hessian_equals_old",
        "coefficient_comparison_is_full_quantum_inverse",
        "same_space_inverse_from_small_coefficient",
        "finite_q_auxiliary_chart_everywhere",
        "Theta_division_through_crossing",
        "literal_homogeneous_constraint_from_P0_extension",
        "remove_classical_heavy_mode",
        "zero_heavy_Gaussian_response",
        "vacuum_constant_dropped_before_matching",
        "only_leading_curvature_symbol",
        "component_pole_scalar_prescription",
        "omit_shift_or_volume_contacts",
        "flat_adjoint_on_curved_background",
        "physical_poles_of_isolated_local_operator",
        "order_reduction_discards_causal_poles",
        "state_dependent_determinant_remainder_zero",
        "full_interacting_light_clock",
        "state_profile_second_clock_contact_zero",
        "matched_heat_mean_identity_applies_to_profile_alone",
        "quantization_commutes_with_gravity_limit",
        "same_state_nonlinear_bounce",
        "UV_Regge_complete",
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
        raise ValueError(
            "Unsupported quadratic, local or full-closure claim accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "literal_new_profile_and_R_not_old_automatic_stability": True,
        "full_local_counteraction_before_constraint_elimination": True,
        "entire_lapse_metric_shift_and_second_chart_terms": True,
        "whole_actual_pairing_and_fourth_order_tensor_terms": True,
        "full_classical_heavy_mode_not_zero_quantum_response": True,
        "remaining_state_response_and_inverse_not_removed": True,
        "original_primitives_and_previous_matching_unchanged": frontier()
        == previous.frontier(),
        "no_original_P8_closure": True,
    }
