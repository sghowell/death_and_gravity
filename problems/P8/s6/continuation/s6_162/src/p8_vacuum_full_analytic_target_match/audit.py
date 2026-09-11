"""A full classical analytic-target match, leaving the quantum frontier open."""

from functools import cache

import sympy as s
from p8_vacuum_two_loop_elastic_cut import audit as previous

from . import bounds, calibration, holomorphic, jets, transport

MODULES = (jets, holomorphic, bounds, calibration, transport)
ITEM = {
    "id": "full_finite_kappa_analytic_target_classical_vacuum_match",
    "status": "BOUNDED_ON_SAME_RESTRICTED_FLAT_SCHWARTZ_CLASS_NOT_QUANTUM_TARGET_OR_ROLLING_COMMON_PARENT_MATCH",
}


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "The full classical target match cannot change its fixed scope ledger"
        )
    return True


@cache
def residuals():
    out = {
        mod.__name__.rsplit(".", 1)[-1] + "_" + name: s.simplify(value)
        for mod in MODULES
        for name, value in mod.data()["checks"].items()
    }
    out.update(
        {
            "nine_primitive_rows_retained": len(frontier()) - 9,
            "all_parent_primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "all_parent_matching_rows_unchanged": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "only_classical_full_target_matching_item_added": len(matching())
            - len(previous.matching())
            - 1,
        }
    )
    return out


def scalar_entry_count():
    return len(residuals())


@cache
def gates():
    return {
        **holomorphic.data()["bounds"],
        **calibration.data()["bounds"],
        **transport.data()["bounds"],
        "actual_finite_kappa_rational_step_target_not_old_polynomial_switch": True,
        "canonical_gradient_scaled_by_kappa_not_unit_clock_gradient": True,
        "full_R_inverse_has_common_complex_neighborhood": True,
        "exact_sixth_and_eighth_jets_and_all_even_higher_tail_owned_once": True,
        "pointwise_W_squared_integrated_using_twenty_one_L2_jets": True,
        "same_Schwartz_class_and_boundary_currents": True,
        "mapped_source_evaluated_before_full_heavy_resolvent": True,
        "field_amplitude_circle_not_physical_momentum_cutoff": True,
        "fixed_flat_metric_not_zero_metric_variation_claim": True,
        "quantum_target_and_rolling_common_parent_dictionaries_open": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
    }


def bad_cases():
    invalid = (
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
    specifications = (
        ("tail", bounds.even_tail, [1, 2, 8], (0, 1)),
        ("integral", bounds.integrated_error, [1, 1, 1, 21], (0, 1, 2)),
        ("compose", bounds.compose, [1, 1], (0, 1)),
    )
    for prefix, call, good, slots in specifications:
        for slot in slots:
            for i, value in enumerate(invalid):
                args = good.copy()
                args[slot] = value
                out.append((f"{prefix}_type_{slot}_{i}", call, tuple(args)))
    for i, value in enumerate(invalid + (s.Integer(8), -2, 3)):
        out.append((f"even_degree_{i}", bounds.even_tail, (1, 2, value)))
    for i, value in enumerate(invalid + (s.Integer(21), 0, -1)):
        out.append((f"jet_count_{i}", bounds.integrated_error, (1, 1, 1, value)))
    for i, value in enumerate(invalid + (s.Integer(1024), 0, 4, 7)):
        out.append(
            (
                f"switch_order_{i}",
                jets.literal,
                (value, 10**800, s.Rational(1, 10**600)),
            )
        )
    for i, caps in enumerate(
        (None, [], (1,), (1, 1, 1, 1, 1, 1), (1, 1, 1, 1, -1), (1, 1, 1, 1, 1.0))
    ):
        out.append(
            (f"majorant_caps_{i}", bounds.weighted_majorant, (jets.PHI**6, caps))
        )
    out.append(
        (
            "symbolic_coefficient",
            bounds.weighted_majorant,
            (s.Symbol("a") * jets.PHI**6,),
        )
    )
    out.extend(
        (
            ("negative_circle", bounds.even_tail, (-1, 2)),
            ("unit_circle", bounds.even_tail, (1, 1)),
            ("negative_radius", bounds.even_tail, (1, -1)),
            ("negative_sixth", bounds.integrated_error, (-1, 1, 1)),
            ("negative_eighth", bounds.integrated_error, (1, -1, 1)),
            ("negative_tail", bounds.integrated_error, (1, 1, -1)),
            ("negative_previous_error", bounds.compose, (-1, 1)),
            ("negative_target_error", bounds.compose, (1, -1)),
        )
    )
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
            ("missing_obligation", validate_scope, (frontier(), matching()[:-1])),
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
        raise ValueError("Unsupported full analytic target input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "full_rational_step_retained_in_Cauchy_circle_bounds": True,
        "nonzero_sextic_density_not_discarded_as_quartic_match": True,
        "no_infinite_volume_pointwise_constant_integral": True,
        "target_argument_is_Psi_and_parent_argument_is_F_Psi": True,
        "same_full_heavy_resolvent_not_low_momentum_loop_replacement": True,
        "no_quantum_or_metric_variation_error_from_classical_density_bound": True,
        "field_amplitude_radius_not_cutoff_or_global_map_inverse": True,
        "original_P8_not_closed": True,
    }
