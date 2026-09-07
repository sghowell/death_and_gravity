"""Separate matrix-ODE, physical-clock and continuous-band Gramian audit.

Synthetic exact solutions test the Taylor mechanism independently of the
physical source model. They are not substituted for its validated replay.
"""

from fractions import Fraction
from functools import cache
from math import factorial

import pytest
import sympy as sp
from flint import arb, arb_mat, ctx
from p8_preparation_cost import core, gramian


@pytest.fixture(autouse=True)
def bounded_arb_context():
    old = ctx.prec, ctx.cap
    ctx.prec, ctx.cap = 256, 24
    try:
        yield
    finally:
        ctx.prec, ctx.cap = old


def exact_arb(expression):
    expression = sp.sympify(expression)
    if expression.is_Rational:
        return arb(int(expression.p))/int(expression.q)
    if expression.is_Add:
        return sum((exact_arb(term) for term in expression.args), arb(0))
    if expression.is_Mul:
        result = arb(1)
        for term in expression.args:
            result *= exact_arb(term)
        return result
    if expression.is_Pow and expression.exp.is_Rational:
        base = exact_arb(expression.base)
        if expression.exp.q == 2:
            return base.sqrt()**int(expression.exp.p)
        if expression.exp.q == 1:
            return base**int(expression.exp.p)
    raise ValueError("The audit only converts exact rational/radical expressions")


@cache
def nilpotent_gramian():
    t, s = sp.symbols("t s", real=True)
    # N(t)=E21+(1+t)E02 has noncommuting constant and linear parts.
    # Its retarded unit-input column is obtained by integrating along the
    # directed chain 1 -> 2 -> 0, not by the production Gramian recurrence.
    first = (t-s)**2/2+(t**3-s**3)/3-s*(t**2-s**2)/2
    impulse = sp.Matrix([first, 1, t-s, 0])
    matrix = (impulse*impulse.T).applyfunc(lambda entry: sp.integrate(entry, (s, 0, t)))
    n0 = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0]])
    n1 = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    return t, matrix, n0, n1


@pytest.mark.parametrize("center", [sp.Rational(0), sp.Rational(1, 3)])
def test_normalized_taylor_recurrence_against_exact_noncommuting_kernel(center):
    t, matrix, n0, n1 = nilpotent_gramian()
    assert n0*n1 != n1*n0
    convert = lambda value: arb_mat(4, 4, [exact_arb(entry) for entry in value])
    generator = [convert(n0+center*n1), convert(n1)]+[arb_mat(4, 4) for _ in range(7)]
    coefficients = gramian._taylor(generator, convert(matrix.subs(t, center)), 8)
    for order, actual in enumerate(coefficients):
        expected = matrix.applyfunc(lambda entry, order=order: sp.diff(entry, t, order).subs(t, center)/factorial(order))
        for i in range(4):
            for j in range(4):
                assert actual[i, j].overlaps(exact_arb(expected[i, j]))
    # Ignoring 1/n! is observably wrong in this exact polynomial control.
    assert sp.diff(matrix[0, 0], t, 3).subs(t, sp.Rational(1, 3)) != sp.diff(matrix[0, 0], t, 3).subs(t, sp.Rational(1, 3))/6


def test_tube_and_lagrange_remainder_enclose_known_exponential_solution():
    # N=I/2, G0=0 implies G11=exp(t)-1. The tube uses lambda=1,
    # ||F||infinity=1 and h=1/4; no guessed solution enclosure is supplied.
    order, h = 3, arb(1)/4
    n0 = arb_mat([[arb(1)/2 if i == j else 0 for j in range(4)] for i in range(4)])
    generator = [n0]+[arb_mat(4, 4) for _ in range(order+1)]
    lam = gramian._infnorm(n0)+gramian._infnorm(n0.transpose())
    tube = h/(1-h*lam)
    assert lam == 1
    assert tube.contains(arb(1)/3)
    range_state = arb_mat(4, 4, [arb(0, tube.upper())]*16)
    remainder = gramian._taylor(generator, range_state, order+1)[-1]
    coefficients = gramian._taylor(generator, arb_mat(4, 4), order)
    polynomial = sum((coefficients[n][1, 1]*h**n for n in range(order+1)), arb(0))
    error = (remainder[1, 1].abs_upper()*h**(order+1)).upper()
    exact = h.exp()-1
    assert not polynomial.overlaps(exact)
    assert (polynomial+arb(0, error)).contains(exact)
    assert tube > exact


@cache
def physical_to_outer_from_kinetic_action():
    u = sp.Symbol("u", negative=True)
    ell, mu = sp.Rational(1, 100), sp.sqrt(39)/2
    d = 1+u**2
    kg, kf = d**6/8, 1/(2*d**6)
    fs, fr = sp.sqrt(2*(kg+kf)), sp.sqrt(2*kg*kf/(kg+kf))
    g, f = sp.Function("g")(u), sp.Function("f")(u)
    light = fs*(kg*g+kf*f)/(kg+kf)
    relative = fr*(f-g)
    fields = sp.Matrix([light, ell*sp.diff(light, u), relative/sp.sqrt(-u),
                        (u*sp.diff(relative, u)-relative/2)/(mu*sp.sqrt(-u))])
    z = sp.symbols("g pg f pf", real=True)
    substitutions = {sp.diff(g, u): z[1]/ell, sp.diff(f, u): z[3]/ell, g: z[0], f: z[2]}
    matrix = fields.subs(substitutions).jacobian(z)
    return matrix.subs(u, -ell).applyfunc(sp.simplify)


def test_endpoint_map_from_literal_kinetic_weights_and_physical_derivatives():
    expected, actual = physical_to_outer_from_kinetic_action(), gramian._endpoint_transform()
    for i in range(4):
        for j in range(4):
            assert actual[i, j].overlaps(exact_arb(expected[i, j]))
    # The incoming radial orientation changes both velocity entries.
    assert expected[3, 1] > 0 and expected[3, 3] < 0
    assert expected[1, 0] != 0 and expected[1, 2] != 0


def test_source_Hilbert_measure_requires_four_ell_cubed_not_fourth_power():
    ell = Fraction(1, 100)
    # In the no-drift calibration X_x=2ell²e2 sigma(a+ell*x), the
    # adjoint loading kernel in du is 2ell e2 on an interval of length ell.
    squared_source_map = (2*ell)**2*ell
    assert squared_source_map == 4*ell**3 == Fraction(4, 10**6)
    assert squared_source_map != 4*ell**4
    assert gramian.certificate()["full_L2_du_scaling"] == squared_source_map


def test_continuous_momentum_bounds_come_from_exact_convex_weights():
    d, w = sp.symbols("d w", positive=True)
    a_k = w/d**4+(1-w)*d**4
    b_k = (1-w)/d**4+w*d**4
    assert sp.factor(d**4-a_k-w*(d**4-d**-4)) == 0
    assert sp.factor(d**4-b_k-(1-w)*(d**4-d**-4)) == 0
    assert sp.expand(sp.Rational(1, 4)-w*(1-w)-(w-sp.Rational(1, 2))**2) == 0
    # d>=1 and 0<w<1 in the physical chart; hence the two diagonal
    # derivatives are <=d^4 and the cross one <=(d^4-d^-4)/2.
    dmax = Fraction(2501, 2500)
    assert dmax**4 < Fraction(101, 100)
    assert (dmax**4-dmax**-4)/2 < Fraction(1, 100)
    cal = core.calibration()
    assert cal["momentum_generator_derivative_entry_sum_upper"] < sp.Rational(1, 5000)
    assert cal["energy_remainder_entry_sum_upper"] < sp.Rational(1, 20)


def test_duration_resolved_Duhamel_and_Frobenius_congruence_constants():
    ell, epsilon_n, epsilon_b = Fraction(1, 100), Fraction(1, 5000), Fraction(1, 2000)
    # The two exponential durations add to1-s, so their product is<2.
    difference = 2*epsilon_n*epsilon_b*Fraction(3, 2)/Fraction(1, 10)
    assert difference == Fraction(3, 10**6)
    # A naive4 instead of2 loses the strict remaining singular margin.
    sigma_mid = Fraction(6, 10**6)
    assert sigma_mid-2*difference == 0
    z_lower = (sigma_mid-difference)**2
    assert z_lower == Fraction(9, 10**12)
    assert z_lower/4 == Fraction(core.calibration()["uniform_Frobenius_gramian_lower"])
    assert ell == Fraction(core.calibration()["ell"])


def test_positive_LDL_does_not_accept_an_indefinite_matrix():
    indefinite = arb_mat([[1, 2, 0, 0], [2, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    with pytest.raises(ValueError, match="positive"):
        gramian._positive_pivots(indefinite)
    pivots = gramian._positive_pivots(arb_mat([[2, 1, 0, 0], [1, 2, 0, 0], [0, 0, 3, 0], [0, 0, 0, 4]]))
    assert all(value > 0 for value in pivots)
