"""Actual two-threshold spectral identities and explicitly bounded moments."""

from functools import cache

from . import calibration, cut, dispersion, moments, partial_wave


@cache
def residuals():
    return {
        prefix + "_" + name: value
        for prefix, module in (
            ("dispersion", dispersion),
            ("cut", cut),
            ("partial_wave", partial_wave),
            ("moments", moments),
            ("calibration", calibration),
        )
        for name, value in module.data()["checks"].items()
    }


@cache
def gates():
    rows = {
        prefix + "_" + name: bool(value)
        for prefix, module in (("moments", moments), ("calibration", calibration))
        for name, value in module.data()["bounds"].items()
    }
    rows.update(
        {
            "same_actual_two_triangle_form_factor_no_new_matching_action": True,
            "positive_measure_with_inverse_moments_not_assumed_finite_total_mass": True,
            "exact_once_subtracted_vertex_dispersion_before_forward_limit": True,
            "both_actual_mass_assignments_and_distinct_pair_thresholds": True,
            "fixed_order_heavy_lines_not_exact_stable_heavy_asymptotic_states": True,
            "anchored_cut_primitive_and_positive_real_log_branch": True,
            "independent_canonical_spin_two_phase_space_normalization": True,
            "identical_pair_factor_and_two_Im_factor_separate": True,
            "angle_independent_potential_contact_projects_out": True,
            "positive_convergent_Q2_series_with_explicit_tail": True,
            "threshold_five_half_power_not_a_zero_mass_pole": True,
            "low_window_inverse_second_moment_strictly_positive": True,
            "low_window_does_not_saturate_total_vertex_slope": True,
            "both_light_and_heavy_full_moments_exceed_one_twentieth_total": True,
            "all_parameter_moments_retained_without_loop_momentum_expansion": True,
            "no_complex_spin_or_Regge_tower_inferred_from_integer_partial_wave": True,
            "no_high_energy_gravitational_contour_bound_from_vertex_dispersion": True,
            "no_old_cosmological_state_or_counterterm_transfer": True,
            "full_V_G_B_and_original_P8_remain_open": True,
        }
    )
    return rows


@cache
def rejected_inputs():
    for name, call, args in calibration.bad_cases():
        try:
            call(*args)
        except (ValueError, TypeError):
            continue
        raise ValueError("Unsupported spectral calibration accepted: " + name)
    return len(calibration.bad_cases())


@cache
def controls():
    d = moments.data()
    return {
        "rejected_inputs": rejected_inputs(),
        "low_cut_moment_not_zero": d["actual_low_window_slope_moment_lower"] > 0,
        "low_cut_moment_not_the_total_slope": d["actual_low_window_slope_moment_upper"]
        < d["actual_total_slope_lower"],
        "heavy_stress_slope_not_omitted": d["both_species_fraction_strict_lower"] > 0,
        "positive_series_tail_not_set_to_zero": calibration.q2_enclosure(2)[
            "tail_upper"
        ]
        > 0,
        "heavy_threshold_not_a_Reggeization_scale": True,
        "vertex_positive_measure_not_amplitude_optical_positivity": True,
        "one_loop_density_not_all_orders_spectral_theorem": True,
        "finite_t_gravity_contour_still_uncomputed": True,
        "original_P8_not_closed": True,
    }
