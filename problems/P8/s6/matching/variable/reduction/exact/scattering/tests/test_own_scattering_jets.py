"""Exact symbolic controls for the independent two-variable Taylor algebra."""

from fractions import Fraction as F

import pytest
import sympy as sp
from p8_exact_stationary.intervals import Interval
from p8_own_scattering.jets import INDICES, TaylorJet, constant, variable


def assert_symbolic_partials(jet, expression, symbols, point):
    for i, j in INDICES:
        exact = sp.diff(expression, symbols[0], i, symbols[1], j).subs(point)
        exact = sp.cancel(exact)
        numerator, denominator = exact.as_numer_denom()
        value = F(int(numerator), int(denominator))
        assert jet.partial(i, j) == Interval(value), (i, j)


@pytest.mark.parametrize("which", ["quotient", "negative_power", "cubic_root"])
def test_all_ten_partials_against_symbolic_differentiation(which):
    v, z = sp.symbols("v z", real=True)
    vj, zj = variable(0, 0), variable(1, 0)
    if which == "quotient":
        expression = (3 + v + 2*z + v*z)**3/(2 + v - z)
        jet = (3 + vj + 2*zj + vj*zj)**3/(2 + vj - zj)
    elif which == "negative_power":
        expression = (-2 + v + 3*z + v*z)**(-3)
        jet = (-2 + vj + 3*zj + vj*zj)**(-3)
    else:
        expression = (8 + v + 2*z + v*z)**sp.Rational(1, 3)
        jet = (8 + vj + 2*zj + vj*zj).cube_root(Interval(2))
    assert_symbolic_partials(jet, expression, (v, z), {v: 0, z: 0})


def test_factorial_normalization_and_total_degree_truncation():
    jet = TaylorJet({(0, 0): 7, (2, 1): F(5, 6)})
    assert jet.coefficient(2, 1) == Interval(F(5, 6))
    assert jet.partial(2, 1) == Interval(F(5, 3))
    v, z = variable(0, 0), variable(1, 0)
    assert all((v**4 + z**4).coefficient(*index) == Interval(0) for index in INDICES)
    assert (v*v*z).partial(2, 1) == Interval(2)


def test_broad_constant_does_not_manufacture_derivative_widths():
    jet = constant(Interval(2, 3)).inverse()
    assert jet.coefficient(0, 0) == Interval(F(1, 3), F(1, 2))
    assert all(jet.coefficient(*index) == Interval(0) for index in INDICES if index != (0, 0))


def test_interval_inverse_derivatives_enclose_entire_positive_base_box():
    jet = variable(0, Interval(2, 3)).inverse()
    assert jet.partial(1, 0) == Interval(F(-1, 4), F(-1, 9))
    assert jet.partial(2, 0) == Interval(F(2, 27), F(1, 4))
    assert jet.partial(3, 0) == Interval(F(-3, 8), F(-2, 27))


def test_exact_inverse_and_cube_root_polynomial_identities_at_a_point():
    v, z = variable(0, 0), variable(1, 0)
    base = 8 + v - z + 2*v*z
    root = base.cube_root(Interval(2))
    inverse_residual = base*base.inverse() - 1
    root_residual = root**3 - base
    assert all(inverse_residual.coefficient(*index) == Interval(0) for index in INDICES)
    assert all(root_residual.coefficient(*index) == Interval(0) for index in INDICES)


@pytest.mark.parametrize("bad", [True, 1.0, float("nan"), float("inf"), "1", None])
def test_inexact_scalar_and_variable_inputs_are_rejected(bad):
    with pytest.raises(TypeError):
        constant(bad)
    with pytest.raises(TypeError):
        variable(0, bad)
    with pytest.raises(TypeError):
        _ = constant(1) + bad


def test_variable_axes_and_indices_are_strict():
    for axis in (True, F(1), 0.0):
        with pytest.raises(TypeError):
            variable(axis, 0)
    with pytest.raises(ValueError):
        variable(2, 0)
    for pair in ((-1, 0), (4, 0), (2, 2)):
        with pytest.raises(ValueError):
            TaylorJet({pair: 1})
        with pytest.raises(ValueError):
            constant(1).partial(*pair)
    with pytest.raises(TypeError):
        TaylorJet({(True, 0): 1})
    with pytest.raises(TypeError):
        TaylorJet({"v": 1})
    with pytest.raises(TypeError):
        TaylorJet([1, 2])


def test_root_and_inverse_chart_guards():
    for base in (0, -8, Interval(-1, 1)):
        with pytest.raises(ValueError):
            constant(base).cube_root(Interval(2))
    for bracket in (Interval(1), Interval(3), Interval(-3, -1)):
        with pytest.raises(ValueError):
            constant(8).cube_root(bracket)
    with pytest.raises(ValueError):
        constant(Interval(7, 9)).cube_root(Interval(2))
    with pytest.raises(TypeError):
        constant(8).cube_root(2.0)
    with pytest.raises(ZeroDivisionError):
        constant(Interval(-1, 1)).inverse()
    for power in (True, F(1, 3), 2.0):
        with pytest.raises(TypeError):
            _ = constant(8)**power


def test_jet_storage_and_export_are_immutable():
    jet = constant(2)
    exported = jet.coefficients()
    exported[(0, 0)] = Interval(9)
    assert jet.coefficient(0, 0) == Interval(2)
    with pytest.raises(AttributeError):
        jet._coefficients = ()
