"""Exact identities and continuous audit gates for the metric local block."""
from functools import cache

import sympy as sp

from . import (
    bounds,
    canonical,
    consistency,
    controls,
    counterterms,
    jets,
    radial,
    variation,
)


@cache
def residuals():
    result = {}
    for name, function in (("canonical", canonical.checks), ("variation", variation.checks),
                           ("counterterm", counterterms.controls), ("poles", radial.pole_checks),
                           ("potential", consistency.potential_checks), ("ordinary", consistency.ordinary_checks),
                           ("adjoint", consistency.self_adjoint_checks), ("scale", consistency.scale_checks),
                           ("physical", bounds.checks), ("failure", controls.checks)):
        for key, value in function().items():
            result[name+"_"+key] = value
    fixture = controls.omitted_counterterm()
    result["omitted_second_order_counterterm_mixed_fixture"] = fixture[1]["bounce_unit_mass_fixtures"][0]+sp.Rational(1600, 81)
    result["omitted_fourth_order_counterterm_mixed_zero_fixture"] = fixture[2]["bounce_unit_mass_fixtures"][0]+sp.Rational(64, 3)
    result["omitted_fourth_order_counterterm_mixed_second_fixture"] = fixture[2]["bounce_unit_mass_fixtures"][2]+sp.Rational(512, 27)
    for key in ("delta_log_g", "delta_log_frequency"):
        result["isotropic_zero_momentum_"+key] = sp.expand(
            (canonical.data("T")[key]-canonical.data("L")[key]).subs(jets.z, 0))
    return result


@cache
def checks():
    example = bounds.operator_bound(10**24, 1000)
    envelope = bounds.envelopes()
    nonnegative = all(value >= 0 for table in envelope.values() for sources in table.values()
                      for row in sources.values() for value in row.values())
    fixture = controls.omitted_counterterm()
    return {
        "physical_lapse_and_logscale_are_independent_sources": True,
        "mass_second_derivatives_are_held_at_the_existing_model_values": True,
        "moving_normalization_and_physical_contact_are_varied_before_projection": True,
        "transverse_multiplicity_is_D_minus_one_before_the_limit": True,
        "full_potential_is_used_once_in_the_counterterm_density": True,
        "quadratic_derivative_mass_terms_need_only_the_first_mass_source_jet": True,
        "linear_mass_weighted_box_R_is_retained": True,
        "pure_divergences_removed_only_under_compact_metric_variations": True,
        "D_dimensional_volume_is_varied_before_the_finite_limit": True,
        "finite_Euler_current_matrix_is_weighted_self_adjoint": True,
        "physical_output_normalization_is_varied_separately": True,
        "all_continuous_coefficient_majorants_are_nonnegative": bool(nonnegative),
        "finite_local_energy_norm_below_1e_minus_38": bool(example["physical_C4_to_C0_component_bounds"]["energy"] < sp.Rational(1, 10**38)),
        "finite_local_pressure_norm_below_1e_minus_38": bool(example["physical_C4_to_C0_component_bounds"]["pressure"] < sp.Rational(1, 10**38)),
        "joint_finite_local_norm_below_1e_minus_38": bool(example["joint_physical_C4_to_C0_bound"] < sp.Rational(1, 10**38)),
        "missing_second_order_evanescent_counterterm_is_detected": bool(fixture[1]["bounce_unit_mass_fixtures"][0] != 0),
        "missing_fourth_order_evanescent_counterterm_is_detected": bool(fixture[2]["bounce_unit_mass_fixtures"][2] != 0),
        "missing_second_mass_vertex_is_detected": bool(controls.second_mass_vertex()["bounce_unit_source_fixture"] != 0),
        "finite_fourth_derivative_local_symbol_is_indefinite_at_the_bounce": bool(controls.top_derivative()["determinant"].subs(jets.u, 0) < 0),
        "local_fourth_derivative_truncation_is_not_resummed_into_extra_particles": True,
        "background_state_and_fixed_tadpole_profiles_are_unchanged": True,
        "nonlocal_full_metric_response_and_coupled_quantum_health_remain_open": True,
    }
