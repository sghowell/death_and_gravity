from fractions import Fraction as Q

import pytest
import sympy as sp
from p8_composite_response import bounds, model
from p8_composite_response.intervals import GRID, Interval, Jet


def test_exact_outward_arithmetic_contains_rational_samples():
    left, right = Interval(Q(-2, 3), Q(7, 5)), Interval(Q(1, 7), Q(3, 2))
    for x in (Q(-2, 3), Q(0), Q(7, 5)):
        for y in (Q(1, 7), Q(1), Q(3, 2)):
            assert (left+right).contains(x+y)
            assert (left*right).contains(x*y)
            assert (left/right).contains(x/y)


def test_radical_integer_enclosure_is_outward():
    for rational in (Q(0), Q(1, 7), Q(2), Q(17, 3), Q(1001, 1000)):
        interval = Interval(rational).sqrt()
        assert interval.lo**2 <= rational <= interval.hi**2
        # The input itself was rounded outward first; below one, sqrt can
        # enlarge that input width before its own outward rounding.
        assert interval.hi-interval.lo <= Q(4, GRID)


def test_zero_denominator_and_inexact_endpoint_controls():
    with pytest.raises(ValueError):
        Interval(-1, 1).reciprocal()
    with pytest.raises(ValueError):
        Interval(-1, 1).sqrt()
    for value in (True, 0.1, float("inf"), float("nan"), sp.Float("0.1")):
        with pytest.raises(TypeError):
            Interval(value)


def test_invalid_power_and_jet_variable_controls():
    for value in (True, Q(1, 2), 1.0):
        with pytest.raises(TypeError):
            Interval(2)**value
        with pytest.raises(TypeError):
            Jet.constant(2)**value
    for value in (-1, 3):
        with pytest.raises(ValueError):
            Jet.variable(1, value)
    with pytest.raises(TypeError):
        Jet.variable(1, True)


def test_first_derivative_jet_against_direct_rational_polynomial():
    x = Jet.variable(Interval(Q(1, 2), Q(3, 2)), 0)
    y = Jet.variable(Interval(2, 3), 1)
    expression = (x*x+3*x*y)/(1+y)
    for xx in (Q(1, 2), Q(1), Q(3, 2)):
        for yy in (Q(2), Q(5, 2), Q(3)):
            assert expression.value.contains((xx*xx+3*xx*yy)/(1+yy))
            assert expression.gradient[0].contains((2*xx+3*yy)/(1+yy))
            assert expression.gradient[1].contains((3*xx-xx*xx)/(1+yy)**2)


def test_model_values_and_derivatives_at_exact_radical_free_points():
    d = model.derive()
    enclosed = bounds.enclosure()
    # R=1 makes X=1 identically at the sample, but derivatives are still
    # taken before substitution. epsilon=0 supplies additional root points.
    points = [(sp.Rational(1, 20000), sp.Rational(1, 3), 2, 1, sp.Rational(1, 100)),
              (0, sp.Rational(2, 3), 4, sp.Rational(3, 2), sp.Rational(1, 128))]
    keys = ("zprime", "Rprime", "kappa", "A_g", "A_f", "V_g", "V_f", "W", "C_lock", "V_lock")
    derivatives = {key: [sp.diff(d[key], variable) for variable in
                         (d["epsilon"], d["z"], d["R"])] for key in keys}
    for eps, u, z, root, q in points:
        at = dict(zip((d["epsilon"], d["u"], d["z"], d["R"], d["q"]),
                      map(sp.sympify, (eps, u, z, root, q)), strict=True))
        for key in keys:
            assert enclosed[key].value.contains(Q(d[key].xreplace(at)))
            for i, derivative in enumerate(derivatives[key]):
                assert enclosed[key].gradient[i].contains(Q(derivative.xreplace(at)))


def test_outward_rounding_is_not_nearest_rounding():
    fraction = Q(1, 3)
    interval = Interval(fraction)
    assert interval.lo < fraction < interval.hi
