from fractions import Fraction as Q
from math import comb

import pytest
import sympy as sp
from p8_variable_global import background


def _multiply(left, right):
    result = [Q(0)]*(len(left)+len(right)-1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i+j] += a*b
    return result


def _fraction_margin_coefficients():
    aa = [Q((comb(10, i) if i <= 10 else 0)
            -(comb(10, i-1) if i else 0), 2) for i in range(12)]
    bb = list(map(Q, (2, 18, 30, 14)))
    ab, b_squared = _multiply(aa, bb), _multiply(bb, bb)
    result = [-value for value in ab]
    for i, value in enumerate(b_squared):
        result[i] += value
    result[0] -= Q(1, 2)
    return result


def test_all_seven_background_equations_with_variable_not_external_clock():
    actual = background.reconstruction_checks()
    assert len(actual) == 7 and not any(actual.values())


def test_fraction_bernstein_reconstruction_compares_every_coefficient():
    coefficients = _fraction_margin_coefficients()
    degree = len(coefficients)-1
    primary = background.bernstein_cells()
    assert degree == 14 and len(primary) == 16
    replay = []
    for cell in range(16):
        power = [sum(coefficients[j]*comb(j, i)*Q(cell)**(j-i)/16**j
                     for j in range(i, degree+1)) for i in range(degree+1)]
        bernstein = [sum(power[i]*Q(comb(j, i), comb(degree, i))
                         for i in range(j+1)) for j in range(degree+1)]
        assert tuple(bernstein) == primary[cell]
        replay.extend(bernstein)
    assert len(replay) == 240 and min(replay) == Q(5, 2)


def test_polynomial_inputs_come_from_the_actual_lapse_comparison():
    u = sp.Symbol("u", real=True)
    d = background.profiles(u)
    assert sp.factor(d["h_prime"]-d["y"]**3*d["A"]) == 0
    bb = sp.diff(d["h"]/d["y"], u)
    assert sp.factor(bb-2*(1+u**2)**2*(1+7*u**2)) == 0
    poly = background.bernstein_margin_polynomial()
    assert sp.factor(bb*(bb-d["A"])-sp.Rational(1, 2)-poly.as_expr().subs(poly.gen, u**2)) == 0


def test_analytic_center_is_not_a_literal_zero_over_zero_or_old_action():
    assert background.center_jets() == {
        "z_prime_0": 1, "z_third_0": 6, "c0": 4, "c_second_0": -16,
        "clock_squared_0": sp.Rational(799, 100),
        "clock_squared_second_0": -sp.Rational(5994, 25)}


def test_clock_bound_uses_the_global_maximum_not_grid_samples():
    u = sp.Symbol("u", real=True)
    aa = background.profiles(u)["A"]
    assert sp.factor(sp.diff(aa, u)-u*(1+u**2)**9*(9-11*u**2)) == 0
    maximum = Q(20, 11)**10/11
    assert maximum < 40 and 1602 < 41**2
    assert Q(1600, 81)-1 == Q(1519, 81) > 0


def test_tail_limits_do_not_keep_a_uniform_positive_algebraic_mass():
    result = background.tail_checks()
    assert result["u22_zprime"] == 1
    assert result["u2_null_stress"] == 8
    assert result["u30_algebraic_mass_squared"] == 4
    z_inf = sp.Symbol("z_infinity", positive=True)
    assert result["u_c"] == 4/z_inf
    assert result["f_Hubble_tail"] == -z_inf


@pytest.mark.parametrize("u", (-3, -1, 0, Q(1, 10), 1, 3))
def test_exact_point_checks_are_only_supplementary_controls(u):
    d = background.profiles(u)
    assert d["z_prime"].is_positive
    assert d["clock_squared"].is_positive
    assert sp.simplify(d["null_stress"]-2*(d["y"]**3*d["z_prime"]-d["h_prime"])) == 0


def test_dropping_the_smoothing_term_loses_positive_clock_at_center():
    aa = Q(1, 2)
    zero_smoothing_zprime = (aa+abs(aa))/2
    assert 2*(8*zero_smoothing_zprime-4)-Q(1, 100) == -Q(1, 100)


@pytest.mark.parametrize("bad", (True, 0.1, sp.Float("0.1"), sp.oo, sp.nan, sp.I))
def test_inexact_or_invalid_profile_inputs_are_not_bounds(bad):
    with pytest.raises((TypeError, ValueError)):
        background.profiles(bad)
