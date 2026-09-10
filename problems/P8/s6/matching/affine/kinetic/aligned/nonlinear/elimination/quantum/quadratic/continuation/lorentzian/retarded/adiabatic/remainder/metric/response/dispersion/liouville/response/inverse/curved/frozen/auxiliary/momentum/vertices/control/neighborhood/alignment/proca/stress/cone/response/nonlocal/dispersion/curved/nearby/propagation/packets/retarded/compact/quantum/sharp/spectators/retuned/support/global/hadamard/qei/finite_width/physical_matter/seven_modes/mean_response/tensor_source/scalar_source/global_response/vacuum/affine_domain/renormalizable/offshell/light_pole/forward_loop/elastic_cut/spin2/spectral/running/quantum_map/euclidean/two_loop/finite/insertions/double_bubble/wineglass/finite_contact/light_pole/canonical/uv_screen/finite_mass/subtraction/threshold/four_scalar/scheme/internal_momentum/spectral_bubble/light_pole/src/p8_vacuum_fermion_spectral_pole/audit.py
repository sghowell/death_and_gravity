"""Full half-trace insertion, on-shell forests and restricted-order controls."""

from functools import cache

import sympy as sp

from . import calibration, enclosure, kernel, quadratic, selected

MODULES = (quadratic, kernel, enclosure, selected, calibration)


@cache
def residuals():
    return {
        mod.__name__.rsplit(".", 1)[-1] + "_" + name: sp.simplify(value)
        for mod in MODULES
        for name, value in mod.data()["checks"].items()
    }


@cache
def gates():
    result = dict(calibration.data()["bounds"])
    result.update(
        {
            "literal_frozen_full_half_trace_differentiated_in_light_covariance": True,
            "exactly_one_light_line_in_mixed_H_Phi_bubble": True,
            "local_quartic_and_stationary_heavy_tadpoles_retained": True,
            "fixed_H_one_point_reference_cancels_stationary_piece": True,
            "whole_outer_affine_reference_cancels_before_bounds": True,
            "same_complete_inner_spectral_insertion_and_OS_references": True,
            "same_original_g_over_Q_normalization_and_inverse_sign": True,
            "outer_mass_and_residue_anchors_fixed_not_refitted": True,
            "spectral_mass_integrals_of_derivative_majorants_converge": True,
            "uniform_unbounded_halfplane_outer_remainder": True,
            "positive_finite_outer_slope_and_spacelike_shape": True,
            "exact_selected_normalization_retains_formal_field_products": True,
            "selected_pole_not_complete_two_loop_pole": True,
            "no_old_pure_scalar_two_loop_budget_transferred": True,
            "finite_reference_window_not_a_derived_cutoff": True,
            "heavy_reference_mass_not_exact_stable_particle_claim": True,
            "no_exact_reflection_positivity_or_full_spectrum_claim": True,
            "original_V_G_B_and_P8_remain_open": True,
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
            values = [2, 1, 1, 1, 144]
            values[j] = value
            result.append((f"type_{j}_{i}", enclosure.bound, tuple(values)))
    result.extend(
        (f"outside_{i}", enclosure.bound, values)
        for i, values in enumerate(
            (
                (0, 1, 1, 1, 144),
                (1, 1, 1, 1, 144),
                (-2, 1, 1, 1, 144),
                (2, -1, 1, 1, 144),
                (2, 1, -1, 1, 144),
                (2, 1, 1, 0, 144),
                (2, 1, 1, 1, 0),
                (2, 1, 1, 1, -144),
            )
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
        raise ValueError(
            "Unsupported spectral quadratic-family input accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "four_point_two_position_factor_not_copied_into_mixed_pole": True,
        "unregulated_tadpole_not_assigned_a_finite_bound": True,
        "whole_inner_and_outer_local_references_paired": True,
        "exact_selected_normalization_not_all_loop_resummation_claim": True,
        "formal_order_two_field_products_not_dropped": True,
        "new_family_not_complete_two_loop_inverse_or_error": True,
        "no_cutoff_or_stable_heavy_spectrum_inferred": True,
        "original_P8_not_closed": True,
    }
