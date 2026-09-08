"""Actual energy weights, derivative ordering and exact integrated tail."""
import pytest
import sympy as sp
from p8_vector_subtraction import controls, tail


def test_exact_formal_subtraction_and_physical_weight_identities():
    assert all(value == 0 for value in tail.checks().values())


def test_both_polarizations_have_the_common_continuous_tail_bound():
    for kind in ("transverse", "longitudinal"):
        d = tail.bounds(kind)
        assert d["reference_subtracted_bracket_over_inverse_frequency_sixth_upper"] < tail.TAIL_CONSTANT


def test_independent_integrated_Proca_control_and_nonzero_fourth_order_UV_term():
    assert all(value == 0 for value in controls.checks().values())
    assert all(value is True for value in controls.proof_checks().values())


def test_UV_divergent_moment_is_not_assigned_a_finite_integral():
    with pytest.raises(ValueError, match="not ultraviolet integrable"):
        controls.finite_moment(sp.Integer(1), 2)


def test_finite_subtracted_integral_and_reference_are_not_full_covariant_matching():
    d = tail.physical_bounds(10**12, 1000)
    assert d["reference_energy_or_pressure_subtracted_tail_over_reference_density"] == sp.Rational(1, 180*10**24)
    assert d["full_subtracted_energy_integral_over_reference_density"] == sp.Rational(4320001, 180*10**24)
    assert d["full_subtracted_pressure_integral_over_reference_density"] == sp.Rational(7776001, 180*10**24)
    assert d["full_covariant_finite_matching_or_all_order_Hadamard_claim"] is False


@pytest.mark.parametrize("bad", (True, 1.0, "1000", sp.oo, 0, 999))
def test_outside_exact_domain_rejected(bad):
    with pytest.raises((TypeError, ValueError)):
        tail.physical_bounds(10**12, bad)
