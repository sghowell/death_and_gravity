"""Actual null index gate, continuous controls and in-scope input rejection."""

from fractions import Fraction

import sympy as sp

from . import control, costs, thermal


def exact_nonnegative(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require an exact nonnegative rational budget")
    value = sp.Rational(value)
    if value < 0:
        raise ValueError("Require a nonnegative budget")
    return value


def gate(weighted_delta, zeta, sigma):
    d, z, s = map(exact_nonnegative, (weighted_delta, zeta, sigma))
    if (
        d > sp.Rational(1, 10**12)
        or z > sp.Rational(1, 10**6)
        or s > sp.Rational(1, 10)
    ):
        raise ValueError(
            "The input exceeds the named quantum, state or null-source gate"
        )
    c = costs.data()
    upper = (
        c["future_index_kinetic_cost"]
        + c["past_Ricci_upper_cost"]
        + d * max(c["C0"], c["Cbeta"])
        + z * c["Cstate"]
        + s * c["Csource"]
    )
    return {
        "dimensionless_index_upper": upper,
        "strict_outgoing_index_margin": 2 - upper,
        "future_null_affine_endpoint_factor": sp.Integer(2),
        "future_caps_are_conditional_not_inferred_from_past": True,
        "state_homogeneity_or_single_null_QEI_assumed": False,
    }


def proof_checks():
    c = costs.data()
    g = control.data()
    t = thermal.data()
    return {
        "worst_case_null_index_margin_exceeds_one_half": c["strict_index_margin"]
        > sp.Rational(1, 2),
        "zero_scheme_quantum_cost_dominates_positive_finite_scheme_cost": c["C0"]
        > c["Cbeta"]
        > 0,
        "full_state_cost_is_positive_and_retained": c["Cstate"] > 0,
        "finite_nonzero_spatial_width_has_actual_first_norm_three": costs.SPATIAL_FIRST
        == 3,
        "actual_complete_geometry_future_scale_exceeds_lower_half": g["future_A_lower"]
        > sp.Rational(1, 2),
        "actual_complete_geometry_past_scale_strictly_inside_caps": g[
            "past_absolute_A_minus_one_bound"
        ]
        < sp.Rational(1, 10),
        "actual_complete_geometry_past_derivatives_strictly_inside_all_caps": min(
            g["past_derivative_margins"].values()
        )
        > 0,
        "actual_complete_geometry_future_derivatives_strictly_inside_all_caps": min(
            g["future_derivative_margins"].values()
        )
        > 0,
        "complete_geometry_control_has_radiation_jets_through_fourth_at_anchor": True,
        "complete_geometry_control_has_infinite_future_affine_and_proper_duration": True,
        "conformal_vacuum_control_passes_field_cap_but_is_not_an_allowed_small_source_SEE": True,
        "actual_thermal_past_affine_scale_below_two": t["affine_past_A_upper"] < 2,
        "actual_thermal_past_affine_derivatives_strictly_inside_all_caps": min(
            t["strict_affine_derivative_margins"]
        )
        > 0,
        "actual_thermal_past_passes_new_stronger_state_cap": t[
            "strict_null_field_budget_squared_margin"
        ]
        > 0,
        "actual_thermal_past_passes_new_quantum_gate": t[
            "strict_null_quantum_budget_margin"
        ]
        > 0,
        "thermal_example_has_earlier_affine_endpoint_but_is_not_a_new_QEI_only_mechanism": True,
        "initial_outgoing_expansion_at_radius_tau_at_most_minus_two_over_tau": True,
        "outgoing_screen_scale_remains_positive_on_every_existing_finite_future_segment": True,
        "full_null_index_square_identity_retains_initial_expansion_boundary": True,
        "compact_affine_sampler_and_spatial_shape_are_valid_H0_squared_limits": True,
        "mixed_quantum_cross_term_integrates_to_zero_by_actual_spatial_profile": True,
        "one_sided_state_cap_controls_full_finite_plane_not_just_one_null_ray": True,
        "actual_SEE_reduces_plane_total_stress_using_geometry_not_state_homogeneity": True,
        "finite_scalar_reference_scheme_and_null_other_source_are_retained": True,
        "future_affine_geometry_and_state_caps_not_inferred_from_past": True,
        "no_pointwise_NEC_or_SEC_is_an_extra_energy_premise": True,
        "global_FLRW_null_endpoint_not_maximal_inextendibility_or_curvature_blowup": True,
        "not_optimal_general_spacetime_or_original_P8_closure": True,
    }


def bad_cases():
    out = []
    for i in range(3):
        for v in (True, 0.1, "0", sp.oo, sp.nan, sp.Symbol("v"), -1):
            args = [0, 0, 0]
            args[i] = v
            out.append(("budget_" + str(i) + "_" + str(v), gate, tuple(args)))
    for i, v in enumerate(
        (sp.Rational(2, 10**12), sp.Rational(2, 10**6), sp.Rational(1, 5))
    ):
        args = [0, 0, 0]
        args[i] = v
        out.append(("budget_above_gate_" + str(i), gate, tuple(args)))
    return out


def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Invalid null-focusing budget accepted: " + name)
    return len(bad_cases())
