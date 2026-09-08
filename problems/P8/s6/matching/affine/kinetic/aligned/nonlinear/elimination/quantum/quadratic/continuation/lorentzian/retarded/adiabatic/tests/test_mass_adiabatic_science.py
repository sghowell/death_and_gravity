"""Varied mode, curved pole, finite-operator and continuous-bound checks."""
import pytest
import sympy as sp
from p8_vector_mass_adiabatic import (
    bounds,
    canonical,
    compact,
    radial,
    variation,
    verify,
)


@pytest.mark.parametrize("function", (canonical.checks, variation.checks, radial.pole_checks, compact.checks, bounds.checks))
def test_exact_science(function):
    verify.affine.certify_residuals(function())


@pytest.mark.parametrize("order", (0, 1, 2))
def test_finite_operator_equals_compact_action(order):
    assert sp.factor(radial.matched(order)-compact.reconstructed(order)) == 0


@pytest.mark.parametrize("order", (1, 2))
def test_missing_finite_counterterm_is_detected(order):
    assert compact.controls()["half_time_unit_mass_adjoint_fixture_"+str(2*order)] != 0


def test_native_order_validation_after_warming_cache():
    radial.laurent(0)
    radial.matched(0)
    compact.coefficients(0)
    for value in (False, 0.0, sp.Integer(0)):
        for function in (radial.laurent, radial.matched, compact.coefficients):
            with pytest.raises(ValueError):
                function(value)


def test_continuous_operator_bound():
    value = bounds.operator_bound(10**24, 1000)["C4_to_C0_upper_bound"]
    assert value == sp.Rational(206506392019063, 236196*10**48)
    assert 0 < value < sp.Rational(1, 10**39)


def test_denominator_validation():
    for value in (1/(1-radial.u**2), 1/(2+radial.u**2), sp.sin(radial.u)):
        with pytest.raises(ValueError):
            bounds.rational_bound(value)


def test_proof_and_rejection_controls():
    assert all(verify.proof_checks().values())
    assert verify.controls()["rejected_inputs"] == 92
