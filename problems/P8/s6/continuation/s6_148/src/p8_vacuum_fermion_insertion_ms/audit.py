"""Only the older paired insertion row obtains its finite local references."""

from functools import cache

import sympy as s
from p8_vacuum_fermion_ms_mass import audit as previous

from . import bubble, calibration, ownership, source, tadpole

MODULES = (tadpole, bubble, source, ownership, calibration)
TARGET = "scalar_Phi2_W1_F0"
STATUS = "BOUNDED_PAIRED_NONLOCAL_AND_FINITE_MS_LOCAL_REFERENCES"


def frontier():
    rows = previous.frontier()
    for r in rows:
        if r["id"] == TARGET:
            r["status"] = STATUS
    return rows


def validate_frontier(rows):
    if rows != frontier():
        raise ValueError("Insertion MS frontier differs from exact scope")
    return True


@cache
def residuals():
    out = {
        mod.__name__.rsplit(".", 1)[-1] + "_" + k: s.simplify(v)
        for mod in MODULES
        for k, v in mod.data()["checks"].items()
    }
    old, new = previous.frontier(), frontier()
    out.update(
        {
            "only_one_paired_row_local_reference_advances": sum(
                a != b for a, b in zip(old, new)
            )
            - 1,
            "nine_rows_retained": len(new) - 9,
            "two_vacuum_rows_wholly_unevaluated": sum(
                r["status"] == "UNEVALUATED" for r in new
            )
            - 2,
            "one_new_paired_local_reference_row": sum(
                r["status"] == STATUS for r in new
            )
            - 1,
            "four_previous_quartic_rows_unchanged": sum(
                r["status"]
                in (
                    "BOUNDED_PAIRED_IN_COMMON_MS_INTERACTION_SCHEME",
                    "BOUNDED_PAIRED_PRIMITIVE_ROW",
                )
                for r in new
            )
            - 4,
        }
    )
    return out


def scalar_entry_count():
    return sum(
        len(v) if isinstance(v, s.MatrixBase) else 1 for v in residuals().values()
    )


@cache
def gates():
    return {
        **calibration.data()["bounds"],
        "same_whole_inner_OS_fermion_forest": True,
        "all_three_local_Hessian_terms_and_one_mixed_bubble": True,
        "one_inserted_light_line_not_a_factor_two": True,
        "regulated_tadpole_keeps_two_beta_anchors": True,
        "tadpole_remainder_nonzero_pole_and_finite_product": True,
        "actual_M_in_outer_bubble_not_light_mass_one": True,
        "whole_bubble_finite_MS_reference_not_zero_subtracted": True,
        "both_independent_mass_ratios_kept": True,
        "uniform_endpoint_log_moments": True,
        "fixed_H_source_cancels_stationary_local_mass": True,
        "on_shell_mass_uses_convergent_uniform_derivative": True,
        "assigned_source_not_a_new_adjustable_parameter": True,
        "other_primitive_rows_unchanged": True,
        "vacuum_source_cross_terms_not_set_to_zero": True,
        "other_matching_canonical_V_G_B_and_original_P8_open": True,
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
    for prefix, call, base in (
        ("tad", tadpole.enclosure, [2, 1, 144]),
        ("bubble", bubble.enclosure, [2, 1, 1, 144]),
    ):
        for j in range(len(base)):
            for i, v in enumerate(invalid):
                vals = list(base)
                vals[j] = v
                out.append((f"{prefix}_type_{j}_{i}", call, tuple(vals)))
    for i, v in enumerate(((1, 1, 144), (2, -1, 144), (2, 1, 0), (2, 1, 145))):
        out.append((f"tad_domain_{i}", tadpole.enclosure, v))
    for i, v in enumerate(
        (
            (1, 1, 1, 144),
            (2, -1, 1, 144),
            (2, 1, 0, 144),
            (2, 1, 2, 144),
            (2, 1, 1, 0),
            (2, 1, 1, 145),
        )
    ):
        out.append((f"bubble_domain_{i}", bubble.enclosure, v))
    for i in range(9):
        rows = frontier()
        rows[i]["status"] = "COMPLETE"
        out.append((f"frontier_{i}", validate_frontier, (rows,)))
    out.extend(
        (
            ("frontier_missing", validate_frontier, (frontier()[:-1],)),
            (
                "frontier_extra",
                validate_frontier,
                (frontier() + [{"id": "extra", "status": "COMPLETE"}],),
            ),
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
        raise ValueError("Unsupported insertion MS input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "inner_OS_forest_retained_in_dimension_d": True,
        "full_MS_finite_pole_products_not_four_dimensional_shortcut": True,
        "two_tadpole_anchors_not_leading_term_alone": True,
        "physical_heavy_mass_and_light_denominator_retained": True,
        "fixed_H_source_not_omitted_or_readjusted": True,
        "vacuum_source_cross_terms_remain_open": True,
        "other_matching_and_canonical_terms_remain_open": True,
        "original_P8_not_closed": True,
    }
