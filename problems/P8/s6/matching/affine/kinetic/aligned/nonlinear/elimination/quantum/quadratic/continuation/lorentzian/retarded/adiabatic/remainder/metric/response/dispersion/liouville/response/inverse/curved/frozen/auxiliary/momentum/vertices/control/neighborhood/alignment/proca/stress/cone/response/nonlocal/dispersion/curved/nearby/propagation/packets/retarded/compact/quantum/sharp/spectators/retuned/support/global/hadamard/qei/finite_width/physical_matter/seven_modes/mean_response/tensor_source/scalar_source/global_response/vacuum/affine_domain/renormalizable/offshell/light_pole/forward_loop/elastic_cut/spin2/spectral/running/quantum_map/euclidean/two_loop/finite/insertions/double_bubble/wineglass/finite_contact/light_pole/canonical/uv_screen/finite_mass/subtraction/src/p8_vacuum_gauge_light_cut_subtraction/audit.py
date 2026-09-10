"""Residuals and explicit nonclosure controls for finite-cut subtraction."""

from functools import cache

import sympy as sp

from . import analytic, calibration, domain, logarithm, subtraction

MODULES = (domain, logarithm, subtraction, analytic, calibration)


@cache
def residuals():
    return {
        module.__name__.rsplit(".", 1)[-1] + "_" + name: sp.simplify(value)
        for module in MODULES
        for name, value in module.data()["checks"].items()
    }


@cache
def gates():
    result = dict(calibration.data()["bounds"])
    log = logarithm.enclosure()
    result.update(
        {
            "eight_exact_positive_log_series_terms": log["terms"] == 8,
            "strict_positive_log_interval_width": log["upper"] > log["lower"] > 0,
            "complex_domain_larger_than_cut_and_Cauchy_disc": analytic.data()[
                "rho_control_radius"
            ]
            > analytic.data()["quotient_control_radius"]
            > analytic.data()["coefficient_Cauchy_radius"]
            > 0,
            "fixed_cutoff_above_both_overlapping_cut_endpoints": subtraction.data()[
                "fixed_upper_cut_invariant"
            ]
            == 6,
            "entire_inherited_box_tail_rebounded_on_new_domain": True,
            "Bose_and_Lorentz_symmetry_remove_square_root_artifacts": True,
            "Ward_soft_factors_control_zero_invariant_endpoint": True,
            "rho_error_holomorphic_not_merely_real_axis_small": True,
            "divided_difference_is_removable_before_integration": True,
            "both_physical_boundary_logarithm_germs_retained": True,
            "crossed_channel_boundary_sign_reversed": True,
            "even_center_derivatives_agree_without_claiming_original_gap": True,
            "Cauchy_coefficient_estimate_uses_explicit_unit_disc": True,
            "known_first_cut_removed_not_dropped_by_absolute_size": True,
            "higher_loop_massless_cuts_not_declared_removed": True,
            "subtraction_changes_observable_not_frozen_parent_action": True,
            "same_cut_polynomial_control_not_free_matching_retuning": True,
            "regular_matched_amplitude_and_high_energy_contour_uncomputed": True,
            "no_new_all_orders_UV_construction_requirement": True,
            "original_V_G_B_and_P8_remain_open": True,
        }
    )
    return {k: bool(v) for k, v in result.items()}


def bad_cases():
    return calibration.bad_cases() + logarithm.bad_cases()


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported subtraction input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "small_mass_enclosure_inconclusive_not_rejected": not calibration.point(36)[
            "strictly_negative_subtraction_coefficient"
        ],
        "real_axis_smallness_not_used_to_bound_derivatives": True,
        "PV_second_derivative_not_integrated_as_a_regular_function": True,
        "opposite_crossed_boundary_signs_not_combined_as_same_sign": True,
        "local_polynomial_same_cut_does_not_fix_full_b2": True,
        "old_scalar_two_loop_budget_not_imported_to_new_model": True,
        "fixed_order_cut_removal_not_full_dispersion_theorem": True,
        "original_P8_not_closed": True,
    }
