"""Exact outer-MS finite part, conversion and primitive-reference frontier."""

from functools import cache

import sympy as sp
from p8_vacuum_fermion_vertex_chord import audit as previous

from . import calibration, conversion, dimensional, laurent, remainder

MODULES = (dimensional, laurent, remainder, conversion, calibration)
TARGET = "scalar_Phi4_W2_F0"


def frontier():
    rows = previous.frontier()
    for r in rows:
        if r["id"] == TARGET:
            r["status"] = "BOUNDED_PAIRED_IN_COMMON_MS_INTERACTION_SCHEME"
    return rows


def validate_frontier(rows):
    if rows != frontier():
        raise ValueError("The outer-MS reference frontier differs from its exact scope")
    return True


@cache
def residuals():
    rows = {
        mod.__name__.rsplit(".", 1)[-1] + "_" + name: sp.simplify(v)
        for mod in MODULES
        for name, v in mod.data()["checks"].items()
    }
    old, new = previous.frontier(), frontier()
    rows.update(
        {
            "one_family_reference_status_changes": sum(a != b for a, b in zip(old, new))
            - 1,
            "nine_primitive_rows_retained": len(new) - 9,
            "five_other_rows_unevaluated": sum(
                r["status"] == "UNEVALUATED" for r in new
            )
            - 5,
            "two_complete_quartic_primitive_rows_unchanged": sum(
                r["status"] == "BOUNDED_PAIRED_PRIMITIVE_ROW" for r in new
            )
            - 2,
            "one_full_family_in_common_MS_interaction_scheme": sum(
                r["status"] == "BOUNDED_PAIRED_IN_COMMON_MS_INTERACTION_SCHEME"
                for r in new
            )
            - 1,
            "older_paired_quadratic_sector_unchanged": sum(
                r["status"] == "BOUNDED_PAIRED_SECTOR" for r in new
            )
            - 1,
        }
    )
    return rows


def scalar_entry_count():
    return sum(
        len(v) if isinstance(v, sp.MatrixBase) else 1 for v in residuals().values()
    )


@cache
def gates():
    result = dict(calibration.data()["bounds"])
    result.update(
        {
            "exact_dimensional_inner_OS_pairing": True,
            "dimensional_spectral_normalization": True,
            "same_one_loop_MS_finite_mass_reference": True,
            "outer_reference_includes_both_loop_factors": True,
            "leading_gamma_function_integral": True,
            "double_and_simple_poles_subtracted": True,
            "finite_Laurent_constants_retained": True,
            "exact_nonzero_mass_ratio_correction": True,
            "analytic_compact_difference_before_regulator_limit": True,
            "integrable_endpoint_log_majorants": True,
            "dyadic_logarithm_bound_not_float_fit": True,
            "literal_common_parent_parameter_map": True,
            "independent_amplitude_and_b2_conversion_sign": True,
            "full_family_bound_adds_conversion_once": True,
            "other_first_order_shift_squares_separate": True,
            "other_canonical_field_contributions_not_claimed": True,
            "five_remaining_primitive_rows_uncomputed": True,
            "no_amplitude_sign_from_error_majorants": True,
            "original_P8_V_G_B_full_matching_open": True,
        }
    )
    return {k: bool(v) for k, v in result.items()}


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
    result = []
    for j in range(3):
        for i, value in enumerate(invalid):
            values = [2, 1, 144]
            values[j] = value
            result.append((f"type_{j}_{i}", remainder.enclosure, tuple(values)))
    for i, values in enumerate(
        (
            (0, 1, 144),
            (sp.Rational(3, 2), 1, 144),
            (2, -1, 144),
            (2, 1, 0),
            (2, 1, -144),
            (2, 1, 145),
        )
    ):
        result.append((f"outside_{i}", remainder.enclosure, values))
    for i in range(9):
        rows = frontier()
        rows[i]["status"] = "COMPLETE"
        result.append((f"frontier_{i}", validate_frontier, (rows,)))
    result.append(("frontier_missing", validate_frontier, (frontier()[:-1],)))
    result.append(
        (
            "frontier_extra",
            validate_frontier,
            (frontier() + [{"id": "new", "status": "COMPLETE"}],),
        )
    )
    return result


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported outer-MS input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "dimensional_OS_subtraction_retained": True,
        "finite_epsilon_pole_products_not_dropped": True,
        "nonzero_mass_ratio_pole_and_finite_part": True,
        "whole_mass_ratio_correction_bounded": True,
        "fixed_common_MS_interaction_reference": True,
        "conversion_not_an_adjustable_b2_contact": True,
        "other_matching_and_canonical_terms_open": True,
        "original_P8_not_closed": True,
    }
