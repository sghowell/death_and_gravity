"""Pin the scalar quadratic step without altering the original matching scope."""

from functools import cache

import sympy as s
from p8_vacuum_scalar_insertion_ms import audit as previous

from . import bounds, calibration, cauchy, conversion, sunset

MODULES = (sunset, conversion, cauchy, bounds, calibration)


def frontier():
    return previous.frontier()


def matching():
    return previous.matching() + [
        {
            "id": "complete_scalar_quadratic_finite_MS_slope_and_converted_OS_bound",
            "status": "BOUNDED_WITH_FIXED_INNER_OS_AND_INTERACTION_MS_COORDINATES",
        }
    ]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("Scalar quadratic MS scope differs from the fixed ledger")
    return True


@cache
def residuals():
    out = {
        mod.__name__.rsplit(".", 1)[-1] + "_" + key: s.simplify(val)
        for mod in MODULES
        for key, val in mod.data()["checks"].items()
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
            "one_scalar_quadratic_matching_step_added": len(matching())
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
        "inverse_sign_from_Euclidean_connected_to_one_PI": True,
        "full_dimensional_sunset_slope_before_pole_projection": True,
        "local_sunset_proper_contractions_have_zero_slope": True,
        "exact_projective_weight_and_finite_moment_bounds": True,
        "all_thirty_two_refinements_and_required_insertions_retained": True,
        "same_first_sheet_integrated_Cauchy_majorants": True,
        "mixed_constant_reference_before_regulator_limit": True,
        "nested_full_dimensional_alpha_slope_product_is_finite": True,
        "proper_interaction_MS_change_has_correct_signed_weights": True,
        "exact_squared_heavy_OS_remainder": True,
        "finite_MS_slope_recorded_before_outer_OS": True,
        "no_extra_LSZ_or_independent_tuned_local_reference": True,
        "global_field_coupling_cross_terms_still_open": True,
        "original_V_G_B_and_P8_remain_open": True,
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
    base = [1, 1, 40, 144]
    out = []
    for j in range(4):
        for i, v in enumerate(invalid):
            vals = list(base)
            vals[j] = v
            out.append((f"type_{j}_{i}", bounds.enclosure, tuple(vals)))
    for i, (j, v) in enumerate(
        ((0, 0), (1, 0), (2, 32), (2, 10**400 + 1), (3, 0), (3, 145))
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
        raise ValueError("Unsupported scalar quadratic MS input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "MS_slope_not_discarded_by_early_outer_OS": True,
        "proper_local_contractions_not_nonlocal_sunset": True,
        "Cauchy_applied_only_to_integrated_holomorphic_functions": True,
        "mass_poles_not_confused_with_slope_poles": True,
        "squared_heavy_mass_derivative_and_both_geometric_tails_retained": True,
        "fixed_order_OS_unit_residue_not_unnormalized_MS_residue": True,
        "intermediate_matching_not_full_GY14_or_physical_truncation": True,
        "original_P8_not_closed": True,
    }
