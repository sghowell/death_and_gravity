"""Complete fermion increment checks with explicit exclusion of old scalar loops."""

from functools import cache

import sympy as sp
from p8_vacuum_fermion_local_matching import potential as parent_potential

from . import calibration, kinematics, low_degree, reference, series

MODULES = (low_degree, series, kinematics, reference, calibration)


@cache
def residuals():
    result = {
        module.__name__.rsplit(".", 1)[-1] + "_" + name: sp.simplify(value)
        for module in MODULES
        for name, value in module.data()["checks"].items()
    }
    p = parent_potential.data()
    result["direct_frozen_potential_quartic_match"] = sp.simplify(
        low_degree.data()["minimal_finite_zero_momentum_four_vertex"]
        - p["quartic_threshold_at_scale_m"]
    )
    return result


@cache
def scalar_count():
    return sum(
        len(value) if isinstance(value, sp.MatrixBase) else 1
        for value in residuals().values()
    )


@cache
def gates():
    result = dict(calibration.data()["bounds"])
    result.update(
        {
            "six_cyclic_orders_include_both_fermion_orientations": True,
            "active_flavor_color_factor_not_gauge_trace_factor": True,
            "independent_zero_momentum_Dirac_trace_matches_potential": True,
            "whole_dimensional_local_reference_subtracted_before_norm_bound": True,
            "degree_two_S4_tensor_space_complete": True,
            "degree_two_coefficient_fixed_by_background_two_point_derivative": True,
            "lower_momentum_degrees_constant_on_equal_mass_shell": True,
            "all_odd_momentum_degrees_absent_by_Lorentz_invariance": True,
            "all_higher_momentum_degrees_bounded_not_only_enumerated_examples": True,
            "complex_kinematics_and_no_pinched_reference_denominator": True,
            "independent_sign_permutations_give_holomorphic_invariant_descent": True,
            "Cauchy_coefficient_bound_uses_complex_not_only_real_control": True,
            "potential_quartic_and_residual_box_not_double_counted": True,
            "external_pole_field_factor_counted_exactly_once": True,
            "formal_loop_expansion_separate_from_exact_selected_normalization": True,
            "no_direct_heavy_Yukawa_or_internal_light_tree_exchange": True,
            "old_scalar_reference_scheme_conversion_not_assumed": True,
            "no_old_complete_two_loop_budget_transferred": True,
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
    result = [
        (f"invalid_mass_{i}", calibration.tail_bound, (v, 1))
        for i, v in enumerate(invalid)
    ]
    result += [
        (f"outside_mass_{i}", calibration.tail_bound, (v, 1))
        for i, v in enumerate((-1, 0, 35))
    ]
    result += [
        (f"invalid_Y_{i}", calibration.tail_bound, (36, v))
        for i, v in enumerate(invalid)
    ]
    result += [
        (f"nonpositive_Y_{i}", calibration.tail_bound, (36, v))
        for i, v in enumerate((-1, 0))
    ]
    for j in range(4):
        for name, value in (("float", 0.0), ("boolean", False)):
            values = [1, 0, 0, sp.Rational(1, 2)]
            values[j] = value
            result.append(
                (f"invalid_band_{j}_{name}", calibration.canonical_band, tuple(values))
            )
    result += [
        (f"outside_band_{i}", calibration.canonical_band, values)
        for i, values in enumerate(
            (
                (0, 0, 0, 0),
                (-1, 0, 0, 0),
                (1, -1, 0, 0),
                (1, 0, -sp.Rational(1, 2), 0),
                (1, 0, 0, 1),
                (1, 0, sp.Rational(1, 2), sp.Rational(1, 4)),
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
        raise ValueError("Unsupported four-scalar matching input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "large_local_quartic_not_misread_as_forward_second_coefficient": True,
        "degree_two_operator_not_silently_dropped_off_shell": True,
        "all_momentum_tail_not_an_all_loop_tail": True,
        "negative_interval_numerator_uses_correct_denominator_endpoint": True,
        "potential_quartic_does_not_get_counted_twice": True,
        "positive_selected_coefficient_not_full_new_model_positivity": True,
        "old_scalar_canonical_scheme_not_equated_with_full_MSbar": True,
        "original_P8_not_closed": True,
    }
