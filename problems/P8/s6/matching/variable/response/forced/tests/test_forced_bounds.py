"""Rational continuous envelopes and declared finite-source error domain."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_variable_forced import bounds, coefficients, exact, independent


@pytest.mark.parametrize("r", [sp.Rational(1, 100), sp.Rational(1, 1000), sp.Rational(1, 10**6)])
def test_separate_fraction_continuous_envelopes(r):
    p = coefficients.envelope(r)
    q = independent.coefficient_constants(Fraction(int(r.p), int(r.q)))
    assert p.keys() == q.keys()
    assert all(p[key] == sp.Rational(q[key]) for key in p)


@pytest.mark.parametrize("r", [sp.Rational(1, 100), sp.Rational(1, 1000), sp.Rational(1, 10**6)])
def test_separate_fraction_causal_norms(r):
    p = bounds.calibration(r)
    q = independent.response_constants(Fraction(int(r.p), int(r.q)))
    assert p.keys() == q.keys()
    assert all(p[key] == sp.Rational(q[key]) for key in p)


def test_monotonicity_is_coefficientwise_not_a_grid():
    checks = independent.checks()
    assert checks["positive_polynomial_coefficients"] > 100
    assert all(value >= 0 for poly in independent.positive_envelope_polynomials().values() for value in poly)


def test_every_norm_margin_strict():
    assert all(value > 0 for value in bounds.checks().values())


def test_physical_source_to_field_error_and_no_delta_convergence_in_bound():
    d1, d2 = sp.Rational(1, 10**12), sp.Rational(1, 10**18)
    assert bounds.physical_response_error(d1, 1, momentum_squared=0) == sp.Rational(1, 10**5)
    assert bounds.physical_response_error(d1, 1) == bounds.physical_response_error(d2, 1)
    assert bounds.physical_response_error(d1, 0) == 0
    assert bounds.physical_response_error(d1, 3) == 3*bounds.physical_response_error(d1, 1)


def test_specific_RMS_proxy_hierarchy_not_cutoff():
    c = bounds.calibration()
    assert c["stiffness_lower_times_r2"] > 9
    assert c["stiffness_endpoint_upper_times_r2"] < sp.Rational(13, 4)**2
    assert c["RMS_over_stiffness_proxy_lower"] == sp.Rational(6, 13)


@pytest.mark.parametrize("bad", [True, 0.01, sp.Float(".01"), 0, -1, sp.Rational(1, 50), sp.oo, sp.nan])
def test_invalid_radius_rejected(bad):
    with pytest.raises((TypeError, ValueError)):
        exact.radius(bad)


@pytest.mark.parametrize("bad", [True, 1e-12, sp.Float("1e-12"), 0, -1, sp.Rational(1, 10**8), sp.oo, sp.nan])
def test_invalid_delta_rejected(bad):
    with pytest.raises((TypeError, ValueError)):
        exact.parameters(bad)


@pytest.mark.parametrize("bad", [True, 0.0, -1, 5, sp.Symbol("K", positive=True), sp.zoo])
def test_invalid_momentum_rejected(bad):
    with pytest.raises((TypeError, ValueError)):
        exact.parameters(sp.Rational(1, 10**12), momentum_squared=bad)


def test_coupled_radius_delta_and_source_domain():
    with pytest.raises(ValueError):
        exact.parameters(sp.Rational(1, 10**9), sp.Rational(1, 10**6))
    for value in (True, -1, 1.0, sp.oo):
        with pytest.raises((TypeError, ValueError)):
            bounds.physical_response_error(sp.Rational(1, 10**12), value)


def test_independent_engine_rejects_inexact_and_negative_radius():
    for value in (True, 0.01, "1/100", -1, 0):
        with pytest.raises((TypeError, ValueError)):
            independent.response_constants(value)
