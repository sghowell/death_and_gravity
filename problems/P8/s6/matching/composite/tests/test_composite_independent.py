from fractions import Fraction

import pytest
from p8_composite import independent as exact


def test_independent_full_coefficient_audit():
    report = exact.checks()
    assert len(report["zero_coefficient_residuals"]) == 11
    assert report["double_root_acceleration"] == "7/36"


def test_direct_chain_rule_on_a_polynomial():
    v = exact.Laurent.variable
    ae, ne = v("Ae"), v("Ne")
    expression = ae**3*ne
    assert not (expression.composite_derivative("a")-3*v("alpha")*ae**2*ne).terms
    assert not (expression.composite_derivative("Nf")-v("beta")*ae**3).terms
    wrong = expression.derivative("a")-3*v("alpha")*ae**2*ne
    assert exact.clear_and_substitute(wrong).terms


def test_composite_substitution_does_not_discard_denominator_numerators():
    v = exact.Laurent.variable
    identity = (v("Ae")-v("alpha")*v("a")-v("beta")*v("b"))/(v("Ae")*v("Ne")**2)
    assert not exact.clear_and_substitute(identity).terms
    assert exact.clear_and_substitute(identity+v("Ne")**-1).terms


def test_direct_fraction_evaluation_and_derivative():
    x = exact.Laurent.variable("a")
    polynomial = Fraction(2, 3)*x**3-2*x+5
    assert polynomial.evaluate({"a": 3}) == 17
    assert polynomial.derivative("a").evaluate({"a": 3}) == 16


@pytest.mark.parametrize("value", [True, 0.5, float("inf"), float("nan"), "1/2"])
def test_inexact_coefficients_rejected(value):
    with pytest.raises(TypeError):
        exact.Laurent.constant(value)


def test_nonmonomial_inverse_rejected():
    with pytest.raises(ValueError, match="monomials"):
        (exact.Laurent.variable("a")+exact.Laurent.variable("b"))**-1
