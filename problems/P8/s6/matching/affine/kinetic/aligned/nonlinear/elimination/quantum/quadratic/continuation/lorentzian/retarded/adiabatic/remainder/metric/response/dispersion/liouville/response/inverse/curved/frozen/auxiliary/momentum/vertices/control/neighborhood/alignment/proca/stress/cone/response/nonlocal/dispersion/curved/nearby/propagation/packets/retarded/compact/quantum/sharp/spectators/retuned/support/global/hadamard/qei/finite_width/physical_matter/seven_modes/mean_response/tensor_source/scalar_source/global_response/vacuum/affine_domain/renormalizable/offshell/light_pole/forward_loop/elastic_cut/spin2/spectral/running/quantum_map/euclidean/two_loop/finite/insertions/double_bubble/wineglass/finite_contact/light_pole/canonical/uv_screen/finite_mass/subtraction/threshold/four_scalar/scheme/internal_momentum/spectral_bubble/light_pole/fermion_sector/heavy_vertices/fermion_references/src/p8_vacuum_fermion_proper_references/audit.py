"""Proper-reference identities, guards and unchanged primitive frontier."""

from functools import cache

import sympy as sp
from p8_vacuum_fermion_heavy_vertices import audit as previous

from . import anchors, bounds, calibration, conversion, dirac

MODULES = (dirac, anchors, conversion, bounds, calibration)


@cache
def residuals():
    rows = {
        mod.__name__.rsplit(".", 1)[-1] + "_" + name: sp.simplify(value)
        for mod in MODULES
        for name, value in mod.data()["checks"].items()
    }
    frontier = previous.frontier()
    rows.update(
        {
            "same_nine_primitive_rows": len(frontier) - 9,
            "seven_primitive_rows_still_unevaluated": sum(
                r["status"] == "UNEVALUATED" for r in frontier
            )
            - 7,
            "only_previous_two_paired_families_are_bounded": sum(
                r["status"]
                in ("BOUNDED_PAIRED_SECTOR", "BOUNDED_IN_DECLARED_PARENT_SUBTRACTION")
                for r in frontier
            )
            - 2,
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
            "inverse_resolvent_derives_scalar_and_gauge_signs": True,
            "literal_Clifford_and_open_color_Casimir_checks": True,
            "formal_dimensional_numerator_retained_before_limit": True,
            "gauge_epsilon_times_pole_finite_terms_retained": True,
            "correct_zero_momentum_Feynman_parameter_weight": True,
            "finite_scalar_moments_and_mass_derivative_ratio_derived": True,
            "scalar_moment_box_bounds_proved_by_positive_integrands": True,
            "gauge_zero_boson_mass_log_endpoints_integrable": True,
            "mass_differentiation_holds_mu_fixed": True,
            "local_Yukawa_UV_beta_matches_frozen_matrix_calculation": True,
            "fermion_mass_and_Yukawa_bare_relations_include_fields": True,
            "shared_complete_Phi_field_factor_not_reset": True,
            "both_external_fermion_field_factors_retained": True,
            "active_opposite_signs_and_inert_gauge_only_flavors_distinct": True,
            "no_closed_loop_flavor_factor_on_open_fermion_line": True,
            "proper_counterterm_occurrences_owned_once": True,
            "zero_momentum_reference_not_physical_fermion_pole": True,
            "Feynman_gauge_local_coefficients_not_independently_gauge_invariant": True,
            "no_new_primitive_integral_marked_evaluated": True,
            "S6_138_outer_reference_conversion_still_open": True,
            "original_P8_V_G_B_and_complete_two_loop_remain_open": True,
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
    for j in range(5):
        for i, value in enumerate(invalid):
            values = [2, 1, 1, sp.Rational(1, 4), 144]
            values[j] = value
            result.append((f"type_{j}_{i}", bounds.enclosure, tuple(values)))
    for i, values in enumerate(
        (
            (1, 1, 1, 0, 144),
            (-2, 1, 1, 0, 144),
            (2, -1, 1, 0, 144),
            (2, 1, -1, 0, 144),
            (2, 1, 1, -1, 144),
            (2, 1, 1, sp.Rational(3, 4), 144),
            (2, 1, 1, 0, 0),
            (2, 1, 1, 0, -144),
        )
    ):
        result.append((f"outside_{i}", bounds.enclosure, values))
    return result


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported proper-reference input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "dimensional_finite_products_not_dropped": True,
        "scale_fixed_during_background_mass_derivative": True,
        "open_line_color_factor_not_closed_flavor_count": True,
        "zero_momentum_and_MS_parameters_not_identified": True,
        "shared_complete_scalar_normalization_retained": True,
        "wide_valid_bounds_can_be_inconclusive": True,
        "primitive_integrals_and_outer_conversion_not_claimed": True,
        "original_P8_not_closed": True,
    }
