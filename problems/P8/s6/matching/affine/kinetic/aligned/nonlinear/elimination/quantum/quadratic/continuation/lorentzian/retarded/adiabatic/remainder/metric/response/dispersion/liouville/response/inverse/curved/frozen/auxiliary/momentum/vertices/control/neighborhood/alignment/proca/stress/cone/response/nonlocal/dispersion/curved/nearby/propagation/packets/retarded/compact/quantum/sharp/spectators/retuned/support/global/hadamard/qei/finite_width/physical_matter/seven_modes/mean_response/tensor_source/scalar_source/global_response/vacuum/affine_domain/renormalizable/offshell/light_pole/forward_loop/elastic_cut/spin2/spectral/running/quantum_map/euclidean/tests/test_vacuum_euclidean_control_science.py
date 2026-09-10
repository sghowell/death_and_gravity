"""Independent all-Euclidean multiplier, restricted-tail and coordinate checks."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_vacuum_euclidean_control import audit, calibration, dyson, kernel, map_domain


@pytest.mark.parametrize(
    "name,value", list(audit.residuals().items()), ids=list(audit.residuals())
)
def test_exact_multiplier_tail_or_coordinate_identity(name, value):
    assert value == 0, name


@pytest.mark.parametrize(
    "name,value", list(audit.gates().items()), ids=list(audit.gates())
)
def test_continuous_domain_or_scope_gate(name, value):
    assert bool(value), name


@pytest.mark.parametrize(
    "name,call,args",
    calibration.bad_cases(),
    ids=[r[0] for r in calibration.bad_cases()],
)
def test_reject_unsupported_Euclidean_control_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "name,value",
    [(k, v) for k, v in audit.controls().items() if k != "rejected_inputs"],
    ids=[k for k in audit.controls() if k != "rejected_inputs"],
)
def test_negative_control_and_original_scope(name, value):
    assert bool(value), name


@pytest.mark.parametrize(
    "y", (0, 1, Fraction(1, 2), sp.Rational(3, 2), 10**400, 10**800)
)
def test_finite_exact_Euclidean_point_enclosure(y):
    d = calibration.point(y)
    assert d["Euclidean_invariant"] == sp.Rational(y)
    assert (
        0
        < d["relative_self_energy_strict_lower"]
        < d["relative_self_energy_strict_upper"]
        < sp.Rational(3, 10**208)
    )
    assert 0 < d["displayed_inverse_relative_lower"] < 1
    assert d["single_line_covariance_relative_upper"] > 1


def test_validation_not_bypassed_by_equal_cached_float():
    calibration.point(1)
    for value in (True, 1.0, sp.Float(1)):
        with pytest.raises((TypeError, ValueError)):
            calibration.point(value)


@pytest.mark.parametrize("order", range(33))
def test_supported_known_kernel_tail_formula(order):
    a = sp.Rational(1, 5)
    d = dyson.tails(a, order)
    assert d["single_line_tail"] == 1 / (1 - a) - sum(a**j for j in range(order + 1))
    assert d["two_line_total_order_tail"] == 1 / (1 - a) ** 2 - sum(
        (j + 1) * a**j for j in range(order + 1)
    )


@pytest.mark.parametrize("order", (0, 1, 2, 4, 8, 16, 32))
def test_actual_insertion_tail_positive_and_decreasing(order):
    d = calibration.insertion_tail(order)
    assert 0 < d["single_line_tail"] < d["two_line_total_order_tail"]
    if order:
        old = calibration.insertion_tail(order - 1)
        assert d["single_line_tail"] < old["single_line_tail"]
        assert d["two_line_total_order_tail"] < old["two_line_total_order_tail"]


def test_unsubtracted_bubble_diverges_not_a_finite_amplitude_norm():
    d = dyson.data()
    assert d["unsubtracted_bubble_diverges_to_positive_infinity"] is True
    assert "not" in d["scope"]


def test_parent_charge_normalization_not_reinterpreted_as_new_gravity_result():
    d = kernel.data()["same_light_and_stress_normalization"]
    assert d["renormalized_stress_charge_not_alpha"] == 1
    assert d["parent_Pi_prime_at_one"] == -d["once_fixed_kinetic_counterterm"]


def test_actual_composite_small_heavy_window_and_large_witness():
    d = map_domain.data()
    assert 0 < d["actual_heavy_window_mixing_upper"] < sp.Rational(1, 10**404)
    assert (
        d["actual_heavy_window_endpoint"]
        < d["explicit_nonuniformity_witness_invariant"]
        < d["selected_Planck_squared"]
    )
    assert d["witness_mixing_strict_rational_lower"] > 1
    assert d["actual_minimum_numerator"] > 0


def test_coordinate_witness_is_prescription_specific_not_physical_cutoff():
    assert "not a physical ghost" in map_domain.data()["scope"]
    assert audit.controls()["large_composite_coordinate_not_a_row_no_go"]


def test_global_known_kernel_tails_not_zero_and_tiny():
    d = calibration.data()
    assert (
        0 < d["actual_single_line_beyond_one_insertion_upper"] < sp.Rational(1, 10**414)
    )
    assert (
        0
        < d["actual_two_line_beyond_one_total_insertion_upper"]
        < sp.Rational(1, 10**414)
    )


def test_exact_counts_and_original_scope():
    assert len(audit.residuals()) == 49
    assert len(audit.gates()) == 34
    assert len(audit.controls()) == 11
    assert audit.rejected_inputs() == 35
