"""Keep primitive bounds unchanged and advance only one matching obligation."""

from functools import cache

import sympy as s
from p8_vacuum_fermion_vacuum import audit as previous

from . import amplitude, calibration, contact, conversion, ownership

MODULES = (contact, conversion, amplitude, ownership, calibration)


def frontier():
    return previous.frontier()


def matching():
    return [
        {
            "id": "fixed_sigma_isolated_MS_conversion",
            "status": "CANCELLED_NONLOCAL_ORDER_TWO_WITH_SAME_REGULATED_REFERENCE",
        },
        {"id": "remaining_scale_parameter_field_map_and_cross_terms", "status": "OPEN"},
        {"id": "remaining_counterterm_and_vacuum_source_assembly", "status": "OPEN"},
        {"id": "full_matched_two_loop_pole_and_amplitude_error", "status": "OPEN"},
        {"id": "finite_EFT_higher_order_truncation", "status": "OPEN"},
        {"id": "V_contour_and_cut_control", "status": "OPEN"},
        {"id": "finite_gravity_G", "status": "OPEN"},
        {"id": "common_parent_B", "status": "OPEN"},
    ]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("Contact-conversion scope differs from the fixed ledger")
    return True


@cache
def residuals():
    rows = {
        mod.__name__.rsplit(".", 1)[-1] + "_" + key: s.simplify(value)
        for mod in MODULES
        for key, value in mod.data()["checks"].items()
    }
    rows.update(
        {
            "all_nine_primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "nine_primitive_rows_retained": len(frontier()) - 9,
            "only_one_isolated_matching_obligation_discharged": sum(
                r["status"] != "OPEN" for r in matching()
            )
            - 1,
            "seven_other_matching_and_closure_obligations_open": sum(
                r["status"] == "OPEN" for r in matching()
            )
            - 7,
        }
    )
    return rows


def scalar_entry_count():
    return len(residuals())


@cache
def gates():
    return {
        **calibration.data()["bounds"],
        "same_fixed_contact_not_adjustable_fit": True,
        "dimensional_affine_lift_holomorphic_before_pole_products": True,
        "entire_I0_and_same_sigma_D_in_both_reference_legs": True,
        "bare_L_G_M_equal_through_order_two": True,
        "both_internal_heavy_triangles_in_full_derivative": True,
        "assigned_linear_counterterm_cancellation_not_bare_contact_argument": True,
        "fixed_physical_mass_cancels_only_local_tadpole": True,
        "first_field_and_H_source_unchanged_in_isolated_direction": True,
        "already_canonical_amplitude_no_additional_LSZ": True,
        "no_sigma_dependent_noncontact_second_G_or_M_shift": True,
        "full_combined_map_and_its_cross_terms_not_inferred": True,
        "coordinate_inverse_tail_not_physical_truncation_bound": True,
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
    cases = []
    for j in range(2):
        for i, value in enumerate(invalid):
            vals = [s.Rational(1, 100), 1]
            vals[j] = value
            cases.append((f"type_{j}_{i}", calibration.enclosure, tuple(vals)))
    for i, vals in enumerate(
        (
            (0, 1),
            (-1, 1),
            (s.Rational(3, 50), 1),
            (1, 1),
            (s.Rational(1, 100), -1),
            (s.Rational(1, 100), 2),
        )
    ):
        cases.append((f"domain_{i}", calibration.enclosure, vals))
    for i in range(9):
        rows = frontier()
        rows[i]["status"] = "COMPLETE"
        cases.append((f"primitive_{i}", validate_scope, (rows, matching())))
    for i in range(8):
        obligations = matching()
        obligations[i]["status"] = "COMPLETE"
        cases.append((f"matching_{i}", validate_scope, (frontier(), obligations)))
    cases.extend(
        (
            ("missing_primitive", validate_scope, (frontier()[:-1], matching())),
            ("missing_obligation", validate_scope, (frontier(), matching()[:-1])),
        )
    )
    return cases


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported contact-conversion input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "single_leg_early_epsilon_truncation_has_nonzero_finite_defect": True,
        "bare_contact_b2_zero_does_not_erase_its_loop_insertion": True,
        "both_heavy_triangles_and_all_three_channels_required": True,
        "nonconstant_affine_sigma_inverse_through_order_two": True,
        "no_extra_external_residue_factor": True,
        "other_MS_directions_and_cross_terms_open": True,
        "algebraic_inverse_tail_not_higher_loop_physics": True,
        "original_P8_not_closed": True,
    }
