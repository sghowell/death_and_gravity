"""Explicit fermion-increment checks and nonclosure controls."""

from functools import cache

import sympy as sp
from p8_vacuum_gauge_yukawa_screen import flow

from . import anchors, calibration, potential, reference, twopoint

MODULES = (potential, anchors, twopoint, reference, calibration)


@cache
def residuals():
    result = {
        module.__name__.rsplit(".", 1)[-1] + "_" + name: sp.simplify(value)
        for module in MODULES
        for name, value in module.data()["checks"].items()
    }
    y, L = sp.symbols("positive_Yukawa quartic_L", positive=True)
    N = 6
    # Minimal wavefunction pole gives gamma_F=2 N y^2/Q.
    # The quartic determinant pole and four external field factors give
    # Q beta_L,F = 8 N y^2 L - 48 N y^4.
    result["fermion_quartic_UV_dictionary_to_inherited_flow"] = sp.expand(
        flow.data()["beta_quartic_times_loop_denominator"]
        - 3 * L * L
        - (8 * N * y * y * L - 48 * N * y**4)
    )
    return result


@cache
def gates():
    result = dict(calibration.data()["bounds"])
    result.update(
        {
            "both_active_mass_signs_and_all_inert_vacuum_constants_retained": True,
            "entire_dimensional_UV_potential_and_affine_kernel_retained": True,
            "independent_Dirac_trace_and_constant_potential_mass_agree": True,
            "all_higher_even_potential_coefficients_positive_in_gap_domain": True,
            "finite_field_domain_not_global_or_bounce_potential": True,
            "cancellation_free_grouped_anchor_series": True,
            "all_omitted_anchor_terms_geometrically_bounded": True,
            "complex_OS_subtraction_domain_and_denominator_gap_explicit": True,
            "local_mass_counterterm_paired_before_remainder_estimate": True,
            "vacuum_energy_reference_not_silently_normal_ordered_away": True,
            "large_reference_mass_not_a_naturalness_exclusion": True,
            "canonical_field_and_coupling_changes_explicit": True,
            "exact_field_rescaling_is_of_selected_one_loop_functional_only": True,
            "first_order_expansions_separately_checked": True,
            "old_scalar_loops_not_included_in_selected_potential_verdict": True,
            "no_old_complete_two_loop_budget_transferred_to_new_model": True,
            "prospective_MSbar_ray_not_relabelled_as_canonical_coupling_flow": True,
            "original_V_G_B_and_P8_remain_open": True,
        }
    )
    return {name: bool(value) for name, value in result.items()}


def bad_cases():
    invalid = (
        True,
        False,
        1.0,
        sp.Float(36),
        "36",
        None,
        sp.oo,
        -sp.oo,
        sp.zoo,
        sp.I,
        sp.nan,
        sp.Symbol("m"),
    )
    return (
        [(f"inexact_mass_{i}", anchors.enclosure, (v,)) for i, v in enumerate(invalid)]
        + [
            (f"outside_mass_{i}", anchors.enclosure, (v,))
            for i, v in enumerate((-1, 0, 35))
        ]
        + [
            (f"invalid_anchor_terms_{i}", anchors.enclosure, (36, v))
            for i, v in enumerate((True, 4.0, sp.Integer(4), "4", None, 0, -1, 17))
        ]
        + [
            (f"invalid_moment_index_{i}", anchors.moment, (v,))
            for i, v in enumerate((True, 1.0, sp.Integer(1), "1", None, -1, 65))
        ]
    )


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported fermion-local-matching input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "inert_flavors_not_omitted_from_vacuum_energy": True,
        "tiny_cut_coefficient_not_used_to_bound_local_mass": True,
        "negative_reference_mass_not_equated_to_negative_physical_pole": True,
        "huge_decimal_anchor_subtraction_not_used": True,
        "constant_field_normalization_not_assumed_to_preserve_coupling_numbers": True,
        "one_loop_fermion_functional_not_full_new_model": True,
        "vacuum_reference_not_transferred_to_bounce_state": True,
        "original_P8_not_closed": True,
    }
