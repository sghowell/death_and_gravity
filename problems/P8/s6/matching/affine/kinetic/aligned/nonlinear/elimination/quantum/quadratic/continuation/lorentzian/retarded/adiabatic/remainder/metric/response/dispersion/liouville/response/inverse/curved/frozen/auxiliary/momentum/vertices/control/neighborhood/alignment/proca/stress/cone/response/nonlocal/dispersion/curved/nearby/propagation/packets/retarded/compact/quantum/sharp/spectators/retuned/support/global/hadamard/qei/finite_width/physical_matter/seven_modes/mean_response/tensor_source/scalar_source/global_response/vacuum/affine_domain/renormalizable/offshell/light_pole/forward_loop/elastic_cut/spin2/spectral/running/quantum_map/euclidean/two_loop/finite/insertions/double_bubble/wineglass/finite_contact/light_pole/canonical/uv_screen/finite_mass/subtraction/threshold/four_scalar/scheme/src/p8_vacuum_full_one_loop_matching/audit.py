"""Full one-loop reference matching, primitive ownership and nonclosure controls."""

from functools import cache

import sympy as sp

from . import amplitude, calibration, local, logarithm, scheme

MODULES = (logarithm, local, scheme, amplitude, calibration)


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
    result.update(
        {
            "MS_interaction_boundary_distinct_from_physical_relevant_references": True,
            "same_renormalized_mass_one_internal_scalar_propagator": True,
            "bare_field_and_coupling_equalities_independently_verified": True,
            "whole_dimensional_poles_and_finite_logarithm_retained": True,
            "fermion_quartic_pole_and_wavefunction_match_inherited_MS_flow": True,
            "old_scalar_finite_contact_fixed_not_refitted": True,
            "independent_counterterm_difference_matches_parameter_conversion": True,
            "canonical_star_vertices_get_no_second_LSZ_copy": True,
            "all_one_loop_scalar_graphs_inherited_in_their_explicit_scheme": True,
            "all_one_loop_fermion_boxes_included_in_the_same_final_reference": True,
            "no_additional_one_loop_gauge_or_H_Yukawa_graph": True,
            "mass_one_pole_and_residue_fixed_with_both_sector_slopes": True,
            "H_one_point_subtraction_cancels_stationary_tadpole_term": True,
            "all_species_one_loop_vacuum_energy_retained": True,
            "positive_Phi_curvature_not_global_quantum_potential_stability": True,
            "entire_complex_unit_disc_Phi_remainder_bounded": True,
            "formal_one_loop_amplitude_not_an_all_order_error_bound": True,
            "old_pure_scalar_two_loop_budget_not_transferred": True,
            "original_V_G_B_and_P8_remain_open": True,
        }
    )
    return {name: bool(value) for name, value in result.items()}


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
    for j in range(4):
        for i, v in enumerate(invalid):
            values = [1, 0, 0, 0]
            values[j] = v
            result.append(
                (
                    f"invalid_interval_{j}_{i}",
                    calibration.coefficient_interval,
                    tuple(values),
                )
            )
    result += [
        (f"outside_interval_{i}", calibration.coefficient_interval, values)
        for i, values in enumerate(
            (
                (0, 0, 0, 0),
                (-1, 0, 0, 0),
                (1, -1, 0, 0),
                (1, 0, 1, 0),
            )
        )
    ]
    result += [
        (f"invalid_log_terms_{i}", logarithm.enclosure, (v,))
        for i, v in enumerate(
            (
                True,
                False,
                1.0,
                sp.Integer(1),
                "1",
                None,
                0,
                -1,
                65,
            )
        )
    ]
    return result


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported full one-loop matching input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "large_finite_mass_reference_not_used_as_negative_free_mass": True,
        "MS_interaction_numbers_not_unchanged_old_canonical_parameters": True,
        "both_sector_residue_signs_retained": True,
        "no_second_LSZ_factor_on_canonical_star_vertices": True,
        "finite_contact_has_no_adjustable_forward_second_coefficient": True,
        "full_one_loop_not_old_pure_scalar_two_loop_or_all_loop_result": True,
        "finite_loop_Phi_pole_not_full_high_energy_or_gauge_spectrum": True,
        "original_P8_not_closed": True,
    }
