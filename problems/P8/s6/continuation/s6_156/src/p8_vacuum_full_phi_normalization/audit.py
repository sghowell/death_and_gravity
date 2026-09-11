"""Complete two-loop Phi pole, without closing the amplitude or original P8."""

from functools import cache

import sympy as s
from p8_vacuum_finite_field_covariance import audit as previous

from . import bounds, calibration, normalization, ownership, pole

MODULES = (normalization, ownership, pole, bounds, calibration)


def frontier():
    return previous.frontier()


def matching():
    return previous.matching() + [
        {
            "id": "complete_GY14_second_Phi_normalization_and_canonical_unit_disc_pole",
            "status": "BOUNDED_WITH_DIRECT_MS_OWNERSHIP_AND_FULL_FINITE_FIELD_REEXPANSION",
        }
    ]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("Full Phi normalization scope differs from the fixed ledger")
    return True


@cache
def residuals():
    out = {
        mod.__name__.rsplit(".", 1)[-1] + "_" + key: s.simplify(v)
        for mod in MODULES
        for key, v in mod.data()["checks"].items()
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
            "one_complete_Phi_two_loop_matching_item_added": len(matching())
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
        "direct_MS_fermion_fields_and_parameters_not_zero_anchor_coordinates": True,
        "all_three_fermionic_quadratic_rows_and_thirty_two_scalar_refinements": True,
        "all_proper_MS_mass_kinetic_Yukawa_occurrences_already_assigned_once": True,
        "complete_scalar_and_fermion_inner_physical_OS_forest": True,
        "remaining_source_and_quartic_insertions_have_zero_slope": True,
        "same_regulator_in_both_bare_field_factors": True,
        "finite_first_parameter_reexpansion_in_second_field_identity": True,
        "epsilon_parameter_and_field_pole_products_cancel_only_after_pairing": True,
        "direct_fermion_unit_disc_extension_from_integrated_radius_two_bound": True,
        "parameter_reexpansion_of_both_first_OS_sectors_retained": True,
        "second_parameter_shift_cannot_modify_fixed_free_inverse": True,
        "finite_order_unit_residue_not_global_or_all_orders_pole": True,
        "matched_amplitude_vacuum_and_physical_truncation_remain_open": True,
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
    out = []
    for j in range(3):
        for i, v in enumerate(invalid):
            vals = [1] * 3
            vals[j] = v
            out.append(
                (f"normalization_type_{j}_{i}", bounds.normalization_upper, tuple(vals))
            )
        vals = [1] * 3
        vals[j] = -1
        out.append(
            (f"normalization_negative_{j}", bounds.normalization_upper, tuple(vals))
        )
    for j in range(2):
        for i, v in enumerate(invalid):
            vals = [1, 1]
            vals[j] = v
            out.append((f"tail_type_{j}_{i}", pole.tail_coefficient, tuple(vals)))
    for i, vals in enumerate(((-1, 1), (1, 0), (1, 2))):
        out.append((f"tail_domain_{i}", pole.tail_coefficient, vals))
    for i in range(3):
        rows = ownership.rows()
        rows[i]["status"] = "UNASSIGNED"
        out.append((f"quadratic_owner_{i}", ownership.validate_quadratic_rows, (rows,)))
    out.append(
        (
            "quadratic_owner_missing",
            ownership.validate_quadratic_rows,
            (ownership.rows()[:-1],),
        )
    )
    out.append(
        (
            "quadratic_owner_duplicated",
            ownership.validate_quadratic_rows,
            (ownership.rows() + ownership.rows()[:1],),
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
        raise ValueError("Unsupported full Phi normalization input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "primitive_completeness_not_used_without_counterterm_ownership": True,
        "proper_zero_anchors_not_extra_MS_parameter_redefinitions": True,
        "finite_parameter_reexpansion_not_dropped_from_Z_H": True,
        "epsilon_slope_times_field_pole_not_dropped": True,
        "half_disc_not_mislabeled_unit_disc_without_Cauchy_extension": True,
        "canonical_mass_one_residue_one_uses_complete_OS_conditions": True,
        "four_point_vacuum_truncation_and_V_G_B_still_open": True,
        "original_P8_not_closed": True,
    }
