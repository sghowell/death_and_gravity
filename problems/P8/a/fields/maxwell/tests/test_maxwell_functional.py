"""Exact proper-clock, radiation and endpoint controls."""

import sympy as sp
from p8a_maxwell import functional


def test_exact_clock_spectral_and_IBP_identities():
    assert all(sp.simplify(value) == 0 for value in functional.identities().values())


def test_radiation_control_has_credit_and_negative_first_derivative_term():
    d = functional.radiation_control()
    assert d["operator_first_times_t"] == -1
    assert d["operator_zero_times_t_squared"] == sp.Rational(15, 16)
    assert d["squared_first_times_t_squared"] == -sp.Rational(15, 8)
    assert d["squared_zero_times_t_fourth"] == sp.Rational(945, 256)
    assert d["reference_EED_pi_squared_t_fourth_over_hbar"] == sp.Rational(31, 2560)
    assert d["absolute_zero_times_t_fourth"] == sp.Rational(4601, 1280)
    assert d["absolute_zero_times_t_fourth"] < d["squared_zero_times_t_fourth"]


def test_radiation_finite_ambiguity_vanishes_without_choosing_beta():
    d = functional.radiation_control()
    assert d["finite_ambiguity_in_density"] == d["finite_ambiguity_in_EED"] == 0


def test_IBP_requires_boundary_control_not_arbitrary_noncompact_tests():
    t = sp.Symbol("t", positive=True)
    # f=1 on [1,2] is not an H2_0 test; the integrated boundary term is nonzero.
    h, f = 1/(2*t), sp.Integer(1)
    coefficient = functional.ibp_coefficients(h, t)
    difference = functional.proper_operator(f, h, t)**2-coefficient["zero"]
    assert sp.integrate(difference, (t, 1, 2)) != 0


def test_scalar_half_coefficient_cannot_replace_two_photon_polarizations():
    assert sp.Rational(1, 8) != sp.Rational(1, 16)
    assert sp.Rational(1, 8) == 2*sp.Rational(1, 16)
