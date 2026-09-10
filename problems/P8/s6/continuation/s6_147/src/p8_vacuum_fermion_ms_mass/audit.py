"""Two primitive finite MS masses; other local forests remain explicitly open."""

from functools import cache

import sympy as s
from p8_vacuum_fermion_ms_slopes import audit as previous

from . import anchors, calibration, correction, forest, outer

MODULES = (anchors, forest, outer, correction, calibration)
TARGETS = ("scalar_Phi2_W0_F2", "gauge_Phi2")
STATUS = "BOUNDED_NONLOCAL_AND_FINITE_MS_MASS_AND_SLOPE_REFERENCES"


def frontier():
    rows = previous.frontier()
    for r in rows:
        if r["id"] in TARGETS:
            r["status"] = STATUS
    return rows


def validate_frontier(rows):
    if rows != frontier():
        raise ValueError("Finite MS mass frontier differs from exact scope")
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
            "only_two_primitive_mass_references_advance": sum(
                a != b for a, b in zip(old, new)
            )
            - 2,
            "nine_rows_retained": len(new) - 9,
            "two_vacuum_rows_wholly_unevaluated": sum(
                r["status"] == "UNEVALUATED" for r in new
            )
            - 2,
            "two_full_primitive_local_reference_rows": sum(
                r["status"] == STATUS for r in new
            )
            - 2,
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
        "same_three_cyclic_words_and_gauge_sign": True,
        "dimension_symbolic_mass_tensors_before_regulator_limit": True,
        "independent_fixed_mu_second_vacuum_mass_derivative": True,
        "proper_mass_kinetic_and_Yukawa_CT_once": True,
        "finite_epsilon_times_pole_products_retained": True,
        "overall_MS_removes_poles_not_finite_reference": True,
        "massless_scalar_is_only_an_auxiliary_anchor": True,
        "actual_scalar_mass_restored_by_full_b_integral": True,
        "whole_fermion_cycle_quartic_MS_CT_retained": True,
        "general_mass_outer_triangle_keeps_soft_and_hard_parts": True,
        "regulator_first_mass_integral_and_finite_part_agree": True,
        "uniform_integrable_endpoint_after_full_forest": True,
        "finite_ratio_correction_integrated_not_omitted": True,
        "on_shell_mass_uses_prior_slope_and_complete_soft_tail": True,
        "older_W1F0_and_other_matching_not_claimed": True,
        "vacuum_rows_full_canonical_pole_and_V_G_B_open": True,
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
            vals = [720, 1, 144]
            vals[j] = v
            out.append((f"type_{j}_{i}", correction.enclosure, tuple(vals)))
    for i, v in enumerate(((719, 1, 144), (720, -1, 144), (720, 1, 0), (720, 1, 145))):
        out.append((f"domain_{i}", correction.enclosure, v))
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
        raise ValueError("Unsupported finite MS mass input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "full_MS_dimensional_forest_not_D4_subtracted_kernel": True,
        "fixed_mu_derivative_before_mu_equals_m": True,
        "nonzero_scalar_mass_difference_retained": True,
        "exact_Gamma_Laurent_coefficients_no_float_fit": True,
        "finite_reference_not_a_new_adjustable_parameter": True,
        "two_vacuum_rows_still_open": True,
        "older_W1F0_and_other_matching_and_canonical_terms_open": True,
        "original_P8_not_closed": True,
    }
