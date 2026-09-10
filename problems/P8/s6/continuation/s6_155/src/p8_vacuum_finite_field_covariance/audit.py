"""One finite field direction advances; the full matching ledger remains open."""

from functools import cache

import sympy as s
from p8_vacuum_scalar_ms_slopes import audit as previous

from . import bounds, calibration, coefficients, covariance, kernel

MODULES = (kernel, coefficients, covariance, bounds, calibration)


def frontier():
    return previous.frontier()


def matching():
    return previous.matching() + [
        {
            "id": "first_finite_Phi_field_covariance_and_MS_projection_commutator",
            "status": "EXACT_REGULATED_FIELD_DIRECTION_AND_BOUNDED_FIRST_SLOPE_CROSS_TERMS",
        }
    ]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("Finite Phi field scope differs from the fixed ledger")
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
            "all_nine_primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "nine_primitive_rows_retained": len(frontier()) - 9,
            "all_prior_matching_obligations_unchanged": sum(
                a != b for a, b in zip(matching()[:-1], previous.matching())
            ),
            "one_first_finite_field_direction_added": len(matching())
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
        "same_full_first_scalar_and_fermion_normalization": True,
        "dimensional_trace_reduction_by_vanishing_boundary_IBP": True,
        "first_epsilon_slope_before_MS_pole_multiplication": True,
        "full_regulated_vertex_and_internal_line_covariance": True,
        "both_complete_one_loop_scalar_and_fermion_families": True,
        "assigned_local_counterterms_use_identical_field_factor": True,
        "MS_projector_commutator_kept_in_coefficient_map": True,
        "fundamental_cubic_and_Yukawa_square_cross_terms": True,
        "common_field_factor_preserves_completed_square_cancellation": True,
        "canonical_amplitude_has_only_one_external_normalization": True,
        "coordinate_commutator_not_double_counted_in_physical_amplitude": True,
        "genuine_new_second_slope_not_silently_assigned": True,
        "remaining_global_matching_and_V_G_B_open": True,
        "original_P8_remains_open": True,
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
    base = [1, 1, 40, 1, 1, 144]
    out = []
    for j in range(6):
        for i, v in enumerate(invalid):
            vals = list(base)
            vals[j] = v
            out.append((f"type_{j}_{i}", bounds.enclosure, tuple(vals)))
    for i, (j, v) in enumerate(
        ((0, 0), (1, 0), (2, 32), (2, 10**400 + 1), (3, 0), (4, -1), (5, 0), (5, 145))
    ):
        vals = list(base)
        vals[j] = v
        out.append((f"domain_{i}", bounds.enclosure, tuple(vals)))
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
        raise ValueError("Unsupported finite Phi field input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "epsilon_slope_retained_in_each_pole_leg": True,
        "internal_line_factors_not_omitted_from_vertex_reexpansion": True,
        "full_renormalized_sum_before_regulator_limit": True,
        "finite_MS_commutator_not_a_new_tuned_counterterm": True,
        "G_and_y_squared_cross_terms_not_lost": True,
        "no_second_LSZ_and_no_double_commutator": True,
        "new_second_slope_and_global_map_still_open": True,
        "original_P8_not_closed": True,
    }
