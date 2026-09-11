"""A clock-transparent coordinate extension, not common-parent closure."""

from functools import cache

import sympy as s
from p8_vacuum_protected_yukawa_profile import audit as previous

from . import calibration, counting, gate, norms, transport

MODULES = (gate, counting, norms, calibration, transport)
ITEM = {
    "id": "polynomial_clock_transparent_physical_source_map_and_full_target_classical_bound",
    "status": "KINEMATIC_CLOCK_ARGUMENT_ALIGNED_WITH_PROTECTED_MASSES_AND_SAME_NAMED_VACUUM_DATA_NOT_BACKGROUND_SOLUTION_QUANTUM_STATE_OR_COMMON_PARENT_CLOSURE",
}


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "The clock-transparent map cannot change the fixed scope ledger"
        )
    return True


@cache
def residuals():
    out = {
        mod.__name__.rsplit(".", 1)[-1] + "_" + key: s.simplify(value)
        for mod in MODULES
        for key, value in mod.data()["checks"].items()
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
            "only_clock_coordinate_matching_item_added": len(matching())
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
        **calibration.data()["bounds"],
        **transport.data()["bounds"],
        "finite_polynomial_gate_no_hidden_infinite_Fourier_support": True,
        "old_and_new_source_supports_distinguished": True,
        "full_heavy_inverse_bounded_before_comparing_actions": True,
        "all_positive_delta_powers_retained_in_linear_majorants": True,
        "common_class_L2_and_Linfinity_norms_not_volume_constants": True,
        "physical_source_action_counterterms_and_Jacobian_transformed_together": True,
        "clock_map_identity_covariant_for_physical_X_equal_one": True,
        "eighth_clock_variation_not_assumed_zero": True,
        "clock_kinematics_not_parent_background_or_quantum_state_solution": True,
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
    for i, value in enumerate(invalid + (-1,)):
        out.append((f"gate_cap_{i}", norms.gate_difference_cap, (value,)))
    for prefix, call in (
        ("exact_action", norms.action_difference),
        ("linear_action", norms.linear_enclosure),
    ):
        good = [s.Rational(1, 10), s.Rational(1, 10), 1, 1, 5000]
        for slot in range(5):
            for i, value in enumerate(invalid):
                args = good.copy()
                args[slot] = value
                out.append((f"{prefix}_type_{slot}_{i}", call, tuple(args)))
        for slot in range(4):
            args = good.copy()
            args[slot] = -1
            out.append((f"{prefix}_negative_{slot}", call, tuple(args)))
        for i, mass in enumerate((0, 4356, 4355)):
            args = good.copy()
            args[-1] = mass
            out.append((f"{prefix}_heavy_gap_{i}", call, tuple(args)))
    out.extend(
        (
            (
                "linear_C_above_quarter",
                norms.linear_enclosure,
                (s.Rational(1, 3), 0, 1, 1, 5000),
            ),
            (
                "linear_delta_above_quarter",
                norms.linear_enclosure,
                (0, s.Rational(1, 3), 1, 1, 5000),
            ),
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
        raise ValueError("Unsupported clock-transparent map input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "finite_gate_not_infinite_analytic_multiplier": True,
        "new_support_thirty_three_and_sixty_six_not_old_three_and_six": True,
        "unexpanded_heavy_resolvent_and_spectral_gap_kept": True,
        "nonzero_new_map_action_error_not_discarded": True,
        "nonzero_eighth_clock_variation_retained": True,
        "clock_identity_not_polynomial_parent_bounce_solution": True,
        "mass_derivative_bounds_not_quantum_state_or_cutoff": True,
        "original_P8_not_closed": True,
    }
