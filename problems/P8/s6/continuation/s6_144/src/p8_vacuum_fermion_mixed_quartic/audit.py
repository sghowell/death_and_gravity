"""Exact mixed-row residuals and mutation-checked nine-row frontier."""

from functools import cache

import sympy as sp
from p8_vacuum_fermion_outer_ms import audit as previous

from . import bounds, catalog, dimensional, joint, reference

MODULES = (catalog, dimensional, joint, reference, bounds)
TARGET = "scalar_Phi4_W1_F2"


def frontier():
    rows = previous.frontier()
    for row in rows:
        if row["id"] == TARGET:
            row["status"] = "BOUNDED_PAIRED_IN_COMMON_MS_INTERACTION_SCHEME"
    return rows


def validate_frontier(rows):
    if rows != frontier():
        raise ValueError("Mixed quartic frontier differs from the exact scope")
    return True


@cache
def residuals():
    result = {
        mod.__name__.rsplit(".", 1)[-1] + "_" + k: sp.simplify(v)
        for mod in MODULES
        for k, v in mod.data()["checks"].items()
    }
    old, new = previous.frontier(), frontier()
    result.update(
        {
            "only_one_row_advances": sum(a != b for a, b in zip(old, new)) - 1,
            "all_nine_rows_retained": len(new) - 9,
            "four_unevaluated_rows_remain": sum(
                r["status"] == "UNEVALUATED" for r in new
            )
            - 4,
            "two_MS_paired_families": sum(
                r["status"] == "BOUNDED_PAIRED_IN_COMMON_MS_INTERACTION_SCHEME"
                for r in new
            )
            - 2,
            "two_chord_primitive_rows_unchanged": sum(
                r["status"] == "BOUNDED_PAIRED_PRIMITIVE_ROW" for r in new
            )
            - 2,
            "one_older_quadratic_paired_sector": sum(
                r["status"] == "BOUNDED_PAIRED_SECTOR" for r in new
            )
            - 1,
        }
    )
    return result


def scalar_entry_count():
    return sum(
        len(v) if isinstance(v, sp.MatrixBase) else 1 for v in residuals().values()
    )


@cache
def gates():
    return {
        **bounds.data()["bounds"],
        "literal_mixed_operator_and_two_assignments": True,
        "all_six_boxes_and_all_heavy_routes": True,
        "proper_quartic_MS_forest_paired_once": True,
        "zero_soft_kernel_is_fixed_mu_mass_derivative": True,
        "signed_spectral_density_not_KL_positivity": True,
        "fermion_only_soft_Cauchy_circle": True,
        "joint_sunset_integrability_at_all_endpoints": True,
        "direct_complex_light_pair_denominator_gap": True,
        "finite_heavy_triangle_moments": True,
        "whole_zero_outer_reference_before_regulator_limit": True,
        "both_soft_and_hard_regulated_reference_pieces": True,
        "finite_pole_products_retained": True,
        "nonzero_mass_ratio_correction_not_discarded": True,
        "literal_finite_parent_parameter_conversion": True,
        "all_four_quartic_primitive_rows_bounded": True,
        "other_matching_conversions_and_quadratic_rows_open": True,
        "V_G_B_and_original_P8_open": True,
    }


def bad_cases():
    invalid = (
        True,
        False,
        1.0,
        sp.Float(1),
        "1",
        None,
        sp.oo,
        -sp.oo,
        sp.zoo,
        sp.I,
        sp.nan,
        sp.Symbol("x"),
    )
    out = []
    for j in range(6):
        for i, value in enumerate(invalid):
            values = [720, 1, 1, 1, 10000, 144]
            values[j] = value
            out.append((f"type_{j}_{i}", bounds.enclosure, tuple(values)))
    for i, values in enumerate(
        (
            (719, 1, 1, 1, 10000, 144),
            (720, -1, 1, 1, 10000, 144),
            (720, 1, -1, 1, 10000, 144),
            (720, 1, 1, -1, 10000, 144),
            (720, 1, 1, 1, 9999, 144),
            (720, 1, 1, 1, 10000, 0),
            (720, 1, 1, 1, 10000, 145),
        )
    ):
        out.append((f"domain_{i}", bounds.enclosure, values))
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
        raise ValueError("Unsupported mixed-row input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "complete_nonlocal_W_retained": True,
        "proper_MS_quartic_and_outer_forest": True,
        "no_large_soft_radius_on_light_lines": True,
        "all_momentum_regions_jointly_integrable": True,
        "finite_epsilon_times_pole_terms_retained": True,
        "not_an_adjustable_forward_contact": True,
        "other_matching_and_canonical_terms_open": True,
        "original_P8_not_closed": True,
    }
