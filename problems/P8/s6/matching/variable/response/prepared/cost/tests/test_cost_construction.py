from fractions import Fraction

import pytest
import sympy as sp
from p8_preparation_cost import construction as c
from p8_preparation_cost import construction_independent as independent


def test_literal_physical_operator_full_boundaries_and_positive_margins():
    result = c.checks()
    assert len(result["residuals"]) == 60
    assert all(value == 0 for value in result["residuals"].values())
    assert all(value > 0 for value in result["strict_margins"].values())


def test_actual_pinned_polarization_source_and_spring_bridge():
    assert len(c.parent_checks()) == 7
    assert set(c.parent_checks().values()) == {0}


def test_complete_independent_function_and_polynomial_replay():
    replay = c.independent_checks()
    assert replay["exact_scalar_function_or_polynomial_comparisons"] == 291
    assert replay["all_equal"]


def test_continuous_tables_and_finite_source_costs():
    values = c.calibration()
    assert values["source_coefficient_bounds"] == (
        (14, 91, 33896, 21745, 10427),
        (3, 4733, 11539, 1139718, 1087204),
        (153, 1258, 606335, 1272779, 56985890))
    assert values["hermite_derivative_bounds"] == (2, 17, 171, 2971, 36005, 499025, 5020519)
    assert values["source_sup_per_scaled_endpoint_linf"] == 445826321
    assert values["source_uu_sup_per_scaled_endpoint_linf"] == 762987669721436
    assert values["source_L2"] == 2*10**8
    assert values["source_second_derivative_L2"] == 4*10**14
    assert values["analytic_gramian_lower"] == sp.Rational(1, 4*10**16)
    assert independent.calibration()["source_L2_raw"] == Fraction(891652642, 5)
    assert not any(isinstance(value, float) for value in independent.calibration().values())


@pytest.mark.parametrize("momentum", [1, sp.Rational(5, 2), 4])
def test_direct_polynomial_endpoint_jet_conditions(momentum):
    d = c.derive()
    endpoint = (sp.Rational(2, 3), -2, sp.Rational(7, 5), sp.Rational(1, 7))
    polynomial = c.hermite_polynomial(endpoint, momentum)
    for order, row in enumerate(c.endpoint_jets()):
        derivative = sp.diff(polynomial, d["x"], order)
        assert derivative.subs(d["x"], 0) == 0
        expected = (row.subs(d["K"], momentum)*sp.Matrix(endpoint))[0]
        assert sp.cancel(derivative.subs(d["x"], 1)-expected) == 0


@pytest.mark.parametrize("matched_order", [3, 4])
def test_missing_source_boundary_jets_break_H02_extension(matched_order):
    # Matching only physical phase data is insufficient for the claimed
    # source regularity. Four jets leave sigma nonzero; five leave sigma'.
    d = c.derive()
    degree = 2*matched_order+1
    coefficients = sp.symbols(f"a0:{matched_order+1}")
    x = d["x"]
    polynomial = sum(value*x**(matched_order+1+j) for j, value in enumerate(coefficients))
    endpoint = sp.Matrix([1, 0, 0, 0])
    jets = [(row.subs(d["K"], 1)*endpoint)[0] for row in c.endpoint_jets()]
    equations = [sp.diff(polynomial, x, j).subs(x, 1)-jets[j] for j in range(matched_order+1)]
    polynomial = polynomial.subs(sp.solve(equations, coefficients))
    assert sp.degree(polynomial, x) <= degree
    derivative_order = matched_order-3
    residual = 0
    for m, value in enumerate(d["source_coefficients"]):
        for j in range(derivative_order+1):
            residual += (sp.binomial(derivative_order, j)
                         *sp.diff(value, d["u"], derivative_order-j).subs({d["u"]: d["right"], d["K"]: 1})
                         *d["ell"]**(-m-j)*sp.diff(polynomial, x, m+j).subs(x, 1))
    assert sp.cancel(residual) != 0


def test_homogeneous_scaling_and_zero_source_are_not_strictly_positive():
    assert all(value == 0 for value in c.source_cost(0).values())
    assert c.source_cost(Fraction(2, 3))["source_L2"] == sp.Rational(4*10**8, 3)
    assert c.calibration()["regular_light_L2_lower"] == 19
    assert c.calibration()["regular_even_L2_lower"] == 900


@pytest.mark.parametrize("invalid", [True, 1.0, sp.Float("0.2"), sp.oo, sp.nan, -1])
def test_source_cost_exact_nonnegative_domain(invalid):
    c.source_cost(1)
    with pytest.raises((TypeError, ValueError)):
        c.source_cost(invalid)


@pytest.mark.parametrize("endpoint,momentum", [
    ((1, 0, 0), 1), ((1.0, 0, 0, 0), 1), ((True, 0, 0, 0), 1),
    ((1, 0, 0, 0), 0), ((1, 0, 0, 0), 5), ((1, 0, 0, 0), 1.0),
])
def test_hermite_exact_endpoint_and_momentum_guards(endpoint, momentum):
    with pytest.raises((TypeError, ValueError)):
        c.hermite_polynomial(endpoint, momentum)


@pytest.mark.parametrize("invalid", [True, False, 0.0])
def test_independent_no_inexact_zero_bypass(invalid):
    with pytest.raises(TypeError):
        independent.Rational({(0, 0): invalid})


def test_independent_polynomial_division_and_jet_guards():
    assert independent.divide_quadratic({(2, 0): Fraction(1), (0, 0): Fraction(1)}, Fraction(1)) == {(0, 0): Fraction(1)}
    assert independent.divide_quadratic({(2, 0): Fraction(1)}, Fraction(1)) is None
    with pytest.raises(ValueError):
        independent.Rational({(0, 0): 1}, True)
    with pytest.raises(ValueError):
        independent.Jet([1])**True
    with pytest.raises(ZeroDivisionError):
        independent.Jet([0, 1]).inverse()
