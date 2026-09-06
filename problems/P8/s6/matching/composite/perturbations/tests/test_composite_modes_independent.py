from fractions import Fraction

import pytest
from p8_composite_modes.independent import Poly, checks, identities, multiply_series


def test_independent_fraction_coefficients_and_normalization_controls():
    result = checks()
    assert len(result["coefficientwise_identities"]) == 8
    assert all(value == "0" for value in result["coefficientwise_identities"].values())
    assert result["printed_4_12_minus_literal_frequency_fixture"] == "-3/2"
    assert result["m_tau_5_mu_u2"] == "1/75"


def test_all_sparse_normal_forms_vanish_without_sampling():
    assert all(not value.terms for value in identities().values())


def test_direct_polynomial_evaluation_and_root_quotient():
    x, y, root = (Poly.variable(name) for name in ("x", "y", "root"))
    expression = (x+y)**3-x**3-y**3-3*x*y*(x+y)
    assert expression.is_zero()
    assert (root**4-Fraction(49, 9)).reduce_root().is_zero()
    assert ((x-2*y)**2).evaluate({"x": Fraction(3, 2), "y": Fraction(1, 4)}) == 1


def test_truncated_exponentials_multiply_to_one():
    plus = [Poly(1), Poly(Fraction(1, 2)), Poly(Fraction(1, 8))]
    minus = [Poly(1), Poly(-Fraction(1, 2)), Poly(Fraction(1, 8))]
    result = multiply_series(plus, minus)
    assert result[0].evaluate({}) == 1
    assert result[1].is_zero() and result[2].is_zero()


@pytest.mark.parametrize("value", [0.1, True, "1/3"])
def test_approximate_coefficient_inputs_rejected(value):
    with pytest.raises(TypeError):
        Poly(value)


def test_omitted_term_is_not_hidden_by_sampled_bounce_zero():
    y = Poly.variable("y")
    omitted = y-1
    assert omitted.evaluate({"y": Fraction(1)}) == 0
    assert not omitted.is_zero()
