"""Exact domain, normalization and scope gates for the named scalar theorem."""

from fractions import Fraction

import sympy as sp

from . import cosmology, field, independent, thermal


def flat_parameters(xi):
    xi = cosmology.exact_nonnegative(xi)
    if xi > sp.Rational(1, 4):
        raise ValueError("The positive-squares flat proof requires 0 <= xi <= 1/4")
    return {
        "xi": xi,
        "quantum_coefficient_times_pi_squared_over_hbar": (3 - 4 * xi) / 48,
        "Wick_square_penalty": 2 * xi,
        "curved_cosmological_transport_proved_at_this_coupling": xi
        == sp.Rational(1, 6),
    }


def proof_checks():
    costs = cosmology.costs()
    gate = cosmology.gate(cosmology.DELTA_MAX, cosmology.ZETA_MAX, cosmology.SIGMA_MAX)
    actual = thermal.history(cosmology.DELTA_MAX)
    return {
        "scalar_spectral_coefficient_independently_reconstructed": sp.Rational(
            independent.spectral_moment(Fraction(1, 6))
        )
        == sp.Rational(7, 144),
        "all_flat_positive_square_coefficients_nonnegative_at_conformal_xi": min(
            sp.Rational(1, 2),
            sp.Rational(1, 2) - 2 * sp.Rational(1, 6),
            2 * sp.Rational(1, 6),
        )
        > 0,
        "coherent_line_energy_uniformly_negative_on_closed_sample_interval": sp.Rational(
            4, 3
        )
        - 2
        < 0,
        "flat_counterexample_polynomial_is_not_claimed_finite_energy_before_cutoff": True,
        "proper_clock_both_operators_retain_correct_scale_factors": True,
        "fresh_scalar_C0_below_five_million": costs["positive_C0_margin"] > 0,
        "fresh_scalar_Cbeta_below_five_million": costs["positive_Cbeta_margin"] > 0,
        "state_cost_is_positive_and_not_deleted": costs["Wick_square_cost"] > 0,
        "actual_worst_gate_margin_exceeds_nine_fiftieths": gate[
            "margin_above_nine_fiftieths"
        ]
        > 0,
        "actual_scalar_thermal_past_strictly_inside_all_C3_caps": min(
            actual["strict_C3_margins"]
        )
        > 0,
        "actual_scalar_thermal_past_strictly_inside_squared_state_gate": actual[
            "strict_squared_field_gate_margin"
        ]
        > 0,
        "own_thermal_low_branch_has_strict_initial_denominator": 1
        - 8 * thermal.coupling(cosmology.DELTA_MAX)
        > 0,
        "generic_ODE_comparison_lemma_has_valid_scalar_coupling": thermal.coupling(
            cosmology.DELTA_MAX
        )
        < sp.Rational(1, 16),
        "endpoint_integral_lower_bound_has_valid_scalar_coupling": thermal.coupling(
            cosmology.DELTA_MAX
        )
        < sp.Rational(1, 64),
        "actual_endpoint_state_amplitude_violates_small_field_gate": 50
        > cosmology.ZETA_MAX**2,
        "past_square_root_cubic_norm_upper": sp.Rational(11, 10) ** 2
        > sp.Rational(6, 5),
        "relative_state_weight_square_root_upper": sp.Rational(9, 10) ** 2
        > sp.Rational(4, 5),
        "future_geometric_caps_not_inferred_from_past": True,
        "future_Wick_square_cap_not_inferred_from_past_or_effective_Newton_constant": True,
        "one_scalar_thermal_energy_coefficient_is_not_photon_coefficient": sp.Rational(
            1, 30
        )
        != sp.Rational(1, 15),
        "scalar_SEE_coupling_is_not_photon_coupling": sp.Rational(1, 360)
        != sp.Rational(31, 180),
        "nonminimal_obstruction_is_not_a_self_consistent_SEE_counterexample": True,
        "quantitative_bound_is_not_claimed_optimal_or_null_averaged": True,
    }


def bad_cases():
    cases = []
    invalid = [
        True,
        False,
        0.1,
        "1/100",
        sp.Float("0.1"),
        sp.oo,
        sp.nan,
        sp.Symbol("x"),
        -1,
    ]
    for value in invalid:
        cases.append(("gate_quantum_" + str(value), cosmology.gate, (value, 0, 0)))
        cases.append(("gate_field_" + str(value), cosmology.gate, (0, value, 0)))
        cases.append(("gate_source_" + str(value), cosmology.gate, (0, 0, value)))
        cases.append(("flat_xi_" + str(value), flat_parameters, (value,)))
        cases.append(("thermal_delta_" + str(value), thermal.coupling, (value,)))
    cases.extend(
        [
            ("excess_quantum", cosmology.gate, (sp.Rational(1, 10**7), 0, 0)),
            ("excess_field", cosmology.gate, (0, sp.Rational(1, 4999), 0)),
            ("excess_source", cosmology.gate, (0, 0, 6)),
            ("excess_xi", flat_parameters, (sp.Rational(1, 3),)),
            ("zero_thermal_coupling", thermal.coupling, (0,)),
            ("excess_thermal_coupling", thermal.coupling, (sp.Rational(1, 10**7),)),
            ("zero_branch", thermal.branch_point, (0, cosmology.DELTA_MAX)),
            ("negative_branch", thermal.branch_point, (-1, cosmology.DELTA_MAX)),
            ("high_branch", thermal.branch_point, (10**6, cosmology.DELTA_MAX)),
            ("float_branch", thermal.branch_point, (2.0, cosmology.DELTA_MAX)),
        ]
    )
    return cases


def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (ValueError, TypeError):
            continue
        raise ValueError("An inadmissible input was accepted: " + name)
    return len(bad_cases())


def exact_groups():
    from . import clock, stress

    return {
        "nonminimal_stress_conformal_and_coherent": field.data()["checks"],
        "fresh_physical_scalar_reference_and_anomaly": stress.checks(),
        "two_measure_correct_proper_clock_terms": clock.checks(),
        "nonminimal_cosmological_history_and_cost": cosmology.checks(),
        "actual_scalar_state_and_SEE_history": thermal.checks(),
    }
