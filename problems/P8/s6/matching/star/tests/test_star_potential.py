from fractions import Fraction

import pytest
import sympy as sp
from p8_star import independent, potential
from p8_star.exact import rational


def test_star_literal_symbolic_identities():
    assert all(value == 0 for value in potential.checks().values())


def test_star_every_coefficient_against_independent_eight_direction_variations():
    for fixture in independent.stress_fixtures():
        actual = potential.evaluate(fixture["beta"], fixture["R"], fixture["N"])
        for key in ("rho_i", "pressure_i", "rho_u", "pressure_u", "null_i", "null_u"):
            assert actual[key] == fixture[key]


def test_star_fraction_polynomial_and_derivative_engine():
    assert all(value.is_zero() for value in independent.polynomial_identities().values())
    report = independent.checks()
    assert report["literal_coframe_jet_directions"] == 24
    assert len(report["actual_algebraic_branch_fixtures"]) == 4


def test_star_lapse_ratio_is_not_the_physical_speed():
    beta = (1, 2, 3, 4, 5)
    correct = potential.evaluate(beta, 2, 3)
    wrong = potential.evaluate(beta, 2, sp.Rational(2, 3))
    assert correct["null_u"] != wrong["null_u"]
    assert correct["null_u"] > 0 and wrong["null_u"] < 0


def test_star_endpoint_terms_change_density_but_not_nulls():
    baseline = potential.evaluate((0, 1, 2, 3, 0), 2, 3)
    endpoints = potential.evaluate((7, 1, 2, 3, -11), 2, 3)
    assert endpoints["rho_i"]-baseline["rho_i"] == 14
    assert endpoints["rho_u"]-baseline["rho_u"] == -22
    assert endpoints["null_i"] == baseline["null_i"]
    assert endpoints["null_u"] == baseline["null_u"]


@pytest.mark.parametrize("value", [True, False, 0.2, sp.Float("0.2"), sp.oo, -sp.oo, sp.nan, sp.zoo,
                                  sp.Symbol("x"), "nan", "1/0"])
def test_star_inexact_nonfinite_or_symbolic_public_inputs_are_rejected(value):
    with pytest.raises((TypeError, ValueError, ZeroDivisionError)):
        rational(value)


def test_star_exact_public_input_forms_are_preserved():
    for value in (Fraction(2, 3), sp.Rational(2, 3), "2/3"):
        assert rational(value) == sp.Rational(2, 3)
    assert rational("0.125") == sp.Rational(1, 8)


@pytest.mark.parametrize("R,N", [(0, 1), (-1, 1), (1, 0), (1, -1)])
def test_star_nonpositive_root_or_lapse_chart_is_rejected(R, N):
    with pytest.raises(ValueError):
        potential.evaluate((0, 1, 2, 3, 4), R, N)


def test_star_wrong_coefficient_count_and_float_beta_are_rejected():
    for beta in ((1, 2), (0, 1, 2, 3, 4, 5), (0, 1.0, 2, 3, 4)):
        with pytest.raises((TypeError, ValueError)):
            potential.evaluate(beta, 1, 1)
