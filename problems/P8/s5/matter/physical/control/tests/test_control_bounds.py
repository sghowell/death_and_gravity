from fractions import Fraction

import pytest
import sympy as sp
from p8_m1_control import bounds, oscillator
from p8_m1_control import model as m


def test_covering_bounds_are_exact_and_nonnegative():
    assert all(Fraction(value) >= 0 for value in bounds.covering_checks().values())


def test_unknown_denominator_and_pole_are_not_certified():
    for value, chart in ((1/(1-2*m.x), "gamma"), (1/m.z, "gamma"),
                         (1/m.LAMBDA, "unitary")):
        with pytest.raises(ValueError, match="Unproved"):
            bounds.coefficient_bound(value, chart)


def test_coefficient_majorant_does_not_sample_its_extremum():
    value = (1-2*m.x+3*m.l*m.z)/(m.J*m.x**2)
    result = bounds.coefficient_bound(value, "unitary")
    assert Fraction(result["absolute_bound"]) == (Fraction(3)+Fraction(3, 10000))*810


def test_uniform_threshold_and_covariant_energy_window():
    result = bounds.build_bounds()
    q = int(result["q_threshold"])
    assert q == 10**20
    assert all(Fraction(value) >= 0 for value in result["threshold_comparison_margins"])
    for row in result["charts"].values():
        assert q >= 2*Fraction(row["mass_operator_bound"])
        assert q >= Fraction(row["covariant_mass_operator_bound"])
    assert Fraction(36, 99) < 1
    assert 1/(1-Fraction(36, 99)) == Fraction(11, 7) < 2
    assert result["not_an_interacting_cutoff_or_optimized_threshold"]


@pytest.mark.parametrize("chart,x_value", [("unitary", sp.Rational(3, 5)),
                                           ("unitary", -sp.Rational(3, 5)), ("gamma", 0)])
def test_high_band_literal_matrices_are_positive(chart, x_value):
    result = oscillator.physical_matrices(chart, x_value, 10**20)
    lower = result["potential"]-sp.Rational(10**20, 2)*sp.eye(2)
    upper = sp.Rational(3*10**20, 2)*sp.eye(2)-result["potential"]
    for matrix in (lower, upper):
        assert matrix[0, 0] > 0
        assert sp.simplify(matrix.det()) > 0
    assert result["connection"].T == -result["connection"]
