from fractions import Fraction

import pytest
from p8_bimetric_general import polynomial as poly
from p8_bimetric_general import verify


def test_coefficientwise_fraction_calculation_and_public_bridge():
    report = poly.exact_checks()
    assert len(report["zero_coefficient_residuals"]) == 13
    assert set(report["zero_coefficient_residuals"].values()) == {0}
    assert all(value == 0 for value in verify.polynomial_bridge().values())


def test_sparse_laurent_derivatives_against_direct_fraction_evaluation():
    a, b = poly.Laurent.variable("a"), poly.Laurent.variable("b")
    expression = (a+b)**3/(2*a**2)
    derivative = expression.derivative("a")
    point = {"a": Fraction(2, 3), "b": Fraction(5, 7)}
    direct = Fraction(3, 2)*(point["a"]+point["b"])**2/point["a"]**2-(point["a"]+point["b"])**3/point["a"]**3
    assert derivative.evaluate(point) == direct
    assert not ((a+b)*(a-b)-a*a+b*b).terms


def test_time_jet_chain_rule_and_omission_controls():
    a, ng = poly.Laurent.variable("a"), poly.Laurent.variable("Ng")
    point = {"a": 2, "ad": 3, "Ng": 5, "Ngd": 7}
    assert (a*a/ng).dt().evaluate(point) == Fraction(12, 5)-Fraction(28, 25)
    assert all(value.terms for value in poly.derive()["controls"].values())


@pytest.mark.parametrize("value", [0.0, 0.1, float("inf"), True])
def test_rounded_or_nonfinite_coefficients_rejected(value):
    with pytest.raises(TypeError, match="exact"):
        poly.Laurent.constant(value)
    with pytest.raises(TypeError, match="exact"):
        poly.Laurent.variable("a")-value


def test_polynomial_inverse_and_zero_lapse_are_not_silently_admitted():
    a, b = poly.Laurent.variable("a"), poly.Laurent.variable("b")
    with pytest.raises(ValueError, match="monomials"):
        (a+b)**-1
    with pytest.raises(ZeroDivisionError):
        (a**-1).evaluate({"a": 0})
