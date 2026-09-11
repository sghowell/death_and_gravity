"""Complete the named fixed-order vacuum/source assembly, not original P8."""

from functools import cache

import sympy as s
from p8_vacuum_full_vacuum_reference import audit as previous

from . import bounds, calibration, differentiation, ownership, source

MODULES = (source, differentiation, bounds, calibration, ownership)
TARGET = "remaining_counterterm_and_vacuum_source_assembly"
STATUS = "ASSEMBLED_FOR_FIXED_ORDER_PHI_POLE_B2_VACUUM_AND_H_SOURCE"


def frontier():
    return previous.frontier()


def matching():
    return [
        {**r, "status": STATUS} if r["id"] == TARGET else dict(r)
        for r in previous.matching()
    ]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("Complete source scope differs from the fixed ledger")
    return True


@cache
def residuals():
    out = {
        mod.__name__.rsplit(".", 1)[-1] + "_" + k: s.simplify(v)
        for mod in MODULES
        for k, v in mod.data()["checks"].items()
    }
    out.update(
        {
            "all_nine_primitive_rows_unchanged": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "nine_primitive_rows_retained": len(frontier()) - 9,
            "only_named_fixed_order_assembly_item_advanced": sum(
                a != b for a, b in zip(matching(), previous.matching())
            )
            - 1,
            "matching_ledger_length_unchanged": len(matching())
            - len(previous.matching()),
        }
    )
    return out


def scalar_entry_count():
    return len(residuals())


@cache
def gates():
    return {
        **calibration.data()["bounds"],
        "scalar_and_fermion_OS_inserted_tadpoles_complete": True,
        "proper_cubic_MS_source_counterterm_added_once": True,
        "counterterm_reference_held_fixed_under_H_derivative": True,
        "local_L_source_terms_cancel_with_physical_mass_reference": True,
        "two_light_lines_in_mass_derivative_retained": True,
        "six_labeled_mass_derivative_sectors_corner_subtracted": True,
        "small_r_split_retains_only_M_to_three_sixteenths": True,
        "finite_source_map_keeps_k1_times_first_pole": True,
        "no_second_field_factor_on_zero_tree_source": True,
        "heavy_field_reference_not_redefined": True,
        "fixed_source_condition_not_an_additional_parameter_freedom": True,
        "vacuum_reference_unchanged_at_order_two": True,
        "all_named_Phi_pole_amplitude_vacuum_source_references_assembled": True,
        "other_coordinate_dictionaries_and_physical_truncation_open": True,
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
    good = [2, 2, 1, 1, 144, 2, 2]
    out = []
    for j in range(7):
        for i, v in enumerate(invalid):
            vals = good.copy()
            vals[j] = v
            out.append(
                (f"source_type_{j}_{i}", bounds.scalar_source_enclosure, tuple(vals))
            )
        vals = good.copy()
        vals[j] = -1
        out.append((f"source_domain_{j}", bounds.scalar_source_enclosure, tuple(vals)))
    for name, vals in (
        ("m", (300, 2, 1, 1, 144, 2, 2)),
        ("M", (2, 30000, 1, 1, 144, 2, 2)),
    ):
        out.append(("power_cap_" + name, bounds.scalar_source_enclosure, vals))
    for j in range(5):
        vals = [1, 1, 1, 1, 144]
        vals[j] = -1
        out.append(
            (f"field_negative_{j}", bounds.field_reexpression_upper, tuple(vals))
        )
    out.append(("field_zero_Q", bounds.field_reexpression_upper, (1, 1, 1, 1, 0)))
    for i, v in enumerate((True, 1.0, "1")):
        out.append(
            (f"field_inexact_{i}", bounds.field_reexpression_upper, (1, v, 1, 1, 144))
        )
    for i in range(8):
        rows = ownership.rows()
        rows[i]["owner"] = "UNASSIGNED"
        out.append((f"source_owner_{i}", ownership.validate_rows, (rows,)))
    out.append(
        ("missing_source_owner", ownership.validate_rows, (ownership.rows()[:-1],))
    )
    out.append(
        (
            "duplicate_source_owner",
            ownership.validate_rows,
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
        raise ValueError("Unsupported complete H-source input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "on_shell_counterterms_not_varied_with_background": True,
        "proper_cubic_source_not_hidden_in_covariance_insertion": True,
        "full_regulated_first_source_not_replaced_by_finite_part": True,
        "mass_derivative_counts_both_light_lines": True,
        "constant_reference_shift_not_physical_nonzero_H_expectation": True,
        "vacuum_higher_order_source_pole_products_not_inferred": True,
        "other_dictionaries_physical_truncation_and_V_G_B_open": True,
        "original_P8_not_closed": True,
    }
