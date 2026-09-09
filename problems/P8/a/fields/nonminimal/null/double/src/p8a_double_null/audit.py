"""Exact sampling budgets, continuous-domain gates and rejected controls."""

from fractions import Fraction

import sympy as sp

from . import boost

DEFAULT_BOOST_RATIO = sp.Rational(4, 3)


def exact(value, *, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require an exact rational, not a float, symbol or limit")
    value = sp.Rational(value)
    if value < 0 or (positive and value == 0):
        raise ValueError("Require a positive width or nonnegative field/coupling")
    return value


def product_budget(
    delta_plus, delta_minus, xi, field_cap_squared, boost_ratio=DEFAULT_BOOST_RATIO
):
    dp = exact(delta_plus, positive=True)
    dm = exact(delta_minus, positive=True)
    r = exact(boost_ratio, positive=True)
    xi, cap = map(exact, (xi, field_cap_squared))
    if xi > sp.Rational(1, 2):
        raise ValueError("The positive-squares proof requires 0<=xi<=1/2")
    h1, h2 = sp.Integer(3), sp.Rational(63, 2)
    cost = (
        sp.Rational(2, 3) * (1 + 4 * xi) * r * h2
        + 4 * (1 - 2 * xi) * h1 * h1 / r
        + sp.Rational(2, 3) * h2 / r**3
    )
    return {
        "quantum_cost_coefficient_of_hbar": sp.factor(
            cost / (8 * sp.pi**2 * dp**3 * dm)
        ),
        "state_cost": sp.factor(8 * xi * cap * h1 / dp**2),
        "boost_squared": r * dp / dm,
        "coupling": xi,
        "single_null_limit_taken": False,
        "curved_or_singularity_bound_assigned": False,
    }


def constants():
    b = boost.data()
    r = sp.Rational(4, 3)
    return {
        "shape_norm_squared": sp.Integer(1),
        "shape_first_norm_squared": b["shape_first_norm_squared"],
        "shape_second_norm_squared": b["shape_second_norm_squared"],
        "actual_shape_derivative_pairings": b[
            "actual_shape_derivative_pairings"
        ].tolist(),
        "actual_product_derivative_Gram": b["actual_product_derivative_Gram"].tolist(),
        "selected_boost_ratio": r,
        "selected_quantum_cost": b["rational_boost_cost"],
        "selected_quantum_coefficient_without_hbar_over_pi_squared": b[
            "rational_boost_cost"
        ]
        / 8,
        "strict_margin_below_74": 74 - b["rational_boost_cost"],
        "selected_state_coefficient": sp.Integer(4),
        "boost_stationary_polynomial_at_1": sp.Integer(35 - 24 - 63),
        "boost_stationary_polynomial_at_4_over_3": 35 * r**4 - 24 * r * r - 63,
    }


def proof_checks():
    c = constants()
    return {
        "chosen_quantum_coefficient_is_strictly_below_74": c["strict_margin_below_74"]
        > 0,
        "unique_positive_boost_stationary_point_lies_between_1_and_4_over_3": c[
            "boost_stationary_polynomial_at_1"
        ]
        < 0
        < c["boost_stationary_polynomial_at_4_over_3"],
        "boost_cost_strictly_convex_and_diverges_at_both_positive_domain_ends": True,
        "positive_squares_coefficients_nonnegative_for_whole_xi_zero_to_half_interval": True,
        "derivative_form_remainder_positive_for_whole_xi_zero_to_half_interval": True,
        "physical_stress_not_square_of_a_transversely_smeared_field": True,
        "Hadamard_pullback_to_timelike_plane_has_spacelike_conormal": True,
        "half_Fourier_domain_keeps_all_spatial_frequencies": True,
        "real_sampler_inversion_keeps_cross_term_in_half_Parseval": True,
        "compact_H0_squared_shape_approached_from_inside_by_smooth_samplers": True,
        "upper_Wick_cap_only_needed_on_full_two_dimensional_sample_region": True,
        "upper_cap_is_not_an_absolute_cap_or_a_momentum_cutoff": True,
        "single_null_shrinking_limit_diverges_and_does_not_contradict_A21": True,
        "curved_transport_only_claimed_for_massless_conformal_scalar": True,
        "curved_vector_is_affinely_parametrized_not_unit_timelike": True,
        "finite_scalar_reference_beta_and_physical_volume_retained": True,
        "no_focusing_theorem_or_SEE_solution_inferred_from_plane_average_alone": True,
        "no_optimal_QEI_or_original_P8_closure_claimed": True,
    }


def bad_cases():
    args = [1, 1, sp.Rational(1, 6), 1, sp.Rational(4, 3)]
    out = []
    for i in range(5):
        for v in (True, 0.5, "1/6", sp.oo, sp.nan, sp.Symbol("v"), -1):
            point = list(args)
            point[i] = v
            out.append(
                ("argument_" + str(i) + "_" + str(v), product_budget, tuple(point))
            )
    for i in (0, 1, 4):
        point = list(args)
        point[i] = 0
        out.append(("zero_positive_argument_" + str(i), product_budget, tuple(point)))
    out.append(("coupling_above_half", product_budget, (1, 1, sp.Rational(3, 4), 1, 1)))
    return out


def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Invalid double-null input accepted: " + name)
    return len(bad_cases())
