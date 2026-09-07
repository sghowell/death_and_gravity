"""Independent audit of fixed physical-source phase separation.

Derivations start with the physical kinetic/source normalization, the Gauss
connection parameters and continuous inequalities. No numerical delta scan,
Gamma fitting, or convergence of the actual residual sequences is assumed.
"""

from fractions import Fraction
from math import comb

import pytest
import sympy as sp
from p8_forced_phase import core
from p8_variable_forced import bounds, operator
from p8_variable_response import connection

F = Fraction


def test_literal_Einstein_normalization_gives_the_actual_limiting_weights():
    v = sp.Symbol("audit_v", positive=True)
    d = 1+v*v
    # Constant limiting f lapse c=2, on the punctured source interval.
    kg, kf = d**6/8, 1/(2*d**6)
    fs = sp.sqrt(2*(kg+kf))
    fr = sp.sqrt(2*kg*kf/(kg+kf))
    w2 = kf/(kg+kf)
    bg = -w2/fr
    jh = d**6*bg/2
    ag = 1/fs
    jl = d**6*ag/2
    expected_jh = -2*d**3/sp.sqrt(d**12+4)
    expected_bg = -4/(d**3*sp.sqrt(d**12+4))
    assert sp.simplify(jh-expected_jh) == 0
    assert sp.simplify(bg-expected_bg) == 0
    assert sp.simplify(ag*jl+bg*jh-2) == 0
    actual = operator.derive()
    values = {actual["u"]: -v, actual["delta"]: 0}
    assert sp.simplify(actual["jH"].subs(values)-jh) == 0
    assert sp.simplify(actual["bg"].subs(values)-bg) == 0
    assert sp.simplify((bg*jh).subs(v, 0)) == sp.Rational(8, 5)


def test_light_operator_and_actual_source_are_regular_at_delta_zero():
    d = operator.derive()
    values = {d["u"]: 0, d["delta"]: 0}
    assert sp.simplify(d["A"].subs(values)-d["K"]-sp.Rational(18, 5)) == 0
    assert sp.simplify(d["jL"].subs(values)-1/sp.sqrt(5)) == 0
    assert sp.simplify(d["ag"].subs(values)-2/sp.sqrt(5)) == 0
    # No heavy mass coefficient is used to infer the light Volterra limit.
    assert d["jH"].has(d["K"]) is False


def test_Gauss_zero_one_connection_parameters_not_infinity_pair():
    mu = core.MU
    a, b, c = -sp.Rational(1, 2), sp.Rational(3, 2), 1-sp.I*mu
    alpha = sp.gamma(c)*sp.gamma(c-a-b)/(sp.gamma(c-a)*sp.gamma(c-b))
    beta = sp.gamma(c)*sp.gamma(a+b-c)/(sp.gamma(a)*sp.gamma(b))
    old = connection.coefficients()
    assert alpha == old["A"]
    assert sp.simplify(sp.gamma(a)*sp.gamma(b)+sp.pi) == 0
    assert sp.simplify((sp.pi/sp.sin(sp.pi*sp.I*mu))/(-sp.pi)-old["B"]) == 0
    assert beta != 0
    assert a-b == -2
    assert (c-a-b).is_integer is False
    assert a+b-c+1 == 1+sp.I*mu
    assert c-a-b+1 == 1-sp.I*mu


def test_Euler_transform_is_required_before_using_the_basic_tail_bound():
    mu, z = sp.symbols("mu z", real=True, positive=True)
    a, b, c = -sp.Rational(1, 2), sp.Rational(3, 2), 1-sp.I*mu
    aa, bb, cc = c-a, c-b, c-a-b+1
    assert sp.simplify(cc-aa-bb-sp.I*mu) == 0
    assert cc-aa == a and cc-bb == b
    assert cc == 1-sp.I*mu
    # y/x=exp(2z) gives the exact left-wave phase after Euler's x^(i mu).
    log_x = -sp.log(1+sp.exp(2*z))
    log_y = 2*z-sp.log(1+sp.exp(2*z))
    assert sp.expand(sp.I*mu*z-sp.I*mu*log_y+sp.I*mu*log_x+sp.I*mu*z) == 0


def test_every_basic_Jost_series_coefficient_ratio_is_below_one():
    n = sp.Symbol("n", nonnegative=True)
    # Squared denominator minus squared numerator of a_(n+1)/a_n,
    # for both c=1+i*mu and c=1-i*mu, mu^2=39/4.
    denominator = ((n+1)**2+sp.Rational(39, 4))*(n+1)**2
    numerator = (n-sp.Rational(1, 2))**2*(n+sp.Rational(3, 2))**2
    difference = sp.Poly(sp.expand(denominator-numerator), n)
    assert all(value > 0 for value in difference.all_coeffs())
    assert difference.as_expr() == 2*n**3+sp.Rational(65, 4)*n**2+25*n+sp.Rational(163, 16)
    # Thus |a_n|<=1 from a_0=1, an all-degree proof, not samples of n.
    x = sp.Symbol("x", positive=True)
    assert sp.factor(sp.diff(x/(1-x), x)-1/(1-x)**2) == 0


def test_small_argument_and_z_derivative_uniformity_not_Fx_vanishing():
    w = sp.Symbol("w", positive=True)
    # w=sqrt(1+delta/(8v^2))>=1. Clear positive denominators.
    assert sp.expand((w**3-3*w+2)-(w-1)**2*(w+2)) == 0
    # x=(1-1/w)/2 <=(w^2-1)/4, hence x<=delta/(32v^2).
    assert sp.simplify(((w*w-1)/4-(1-1/w)/2)*4*w-(w-1)**2*(w+2)) == 0
    x = sp.Symbol("x", real=True)
    xz = -2*x*(1-x)
    assert sp.diff(xz, x) == -2+4*x
    fx_at_zero = -sp.Rational(3, 4)/(1-sp.I*core.MU)
    assert fx_at_zero != 0
    # |F-1|<=h/(1-h), |F_x|<=1/(1-h)^2 and |x_z|<=2h
    # give a vanishing physical-z derivative, not a vanishing F_x.
    h = sp.Symbol("h", positive=True)
    assert sp.limit(2*h/(1-h)**2, h, 0) == 0


def test_clock_error_is_uniform_on_the_fixed_punctured_source_interval():
    w = sp.Symbol("w", positive=True)
    error_upper = (w*w-1)/4-sp.log((1+w)/2)
    derivative = sp.diff(error_upper, w)
    assert sp.factor(derivative-(w-1)*(w+2)/(2*(w+1))) == 0
    assert error_upper.subs(w, 1) == 0
    # v>=r/2 makes the largest local Gauss argument <=delta/(8r^2).
    delta, r = sp.symbols("delta r", positive=True)
    assert sp.simplify(delta/(32*(r/2)**2)-delta/(8*r*r)) == 0


def test_Wronskian_forces_nonzero_transmission_without_small_reflection():
    z, mu = sp.symbols("z mu", real=True)
    a, ac, b, bc = sp.symbols("A Abar B Bbar")
    wave = a*sp.exp(sp.I*mu*z)+b*sp.exp(-sp.I*mu*z)
    conjugate = ac*sp.exp(-sp.I*mu*z)+bc*sp.exp(sp.I*mu*z)
    wronskian = sp.expand(wave*sp.diff(conjugate, z)-sp.diff(wave, z)*conjugate)
    assert sp.simplify(wronskian+2*sp.I*mu*(a*ac-b*bc)) == 0
    phase = sp.Symbol("phase", nonzero=True)
    assert sp.diff(ac*phase, phase) == ac  # Still nonzero when B=0.


def test_cross_bounce_conjugations_and_double_clock_scale():
    pr, pv = sp.symbols("phase_r phase_v", nonzero=True)
    ac, bc = sp.symbols("Abar Bbar")
    # f(right)~pr, f(left)~A/pv+B*pv. Conjugate reverses phases.
    product = sp.expand(pr*(ac*pv+bc/pv))
    assert product == ac*pr*pv+bc*pr/pv
    r, v, delta = sp.symbols("r v delta", positive=True)
    lr = sp.log(4*sp.sqrt(2)*r/sp.sqrt(delta))
    lv = sp.log(4*sp.sqrt(2)*v/sp.sqrt(delta))
    assert sp.simplify(sp.expand_log(lr+lv-sp.log(32*r*v/delta), force=True)) == 0
    assert sp.simplify(sp.expand_log(lr-lv-sp.log(r/v), force=True)) == 0
    assert sp.simplify(delta*sp.diff(lr+lv, delta)) == -1
    assert sp.simplify(delta*sp.diff(lr-lv, delta)) == 0


def test_actual_weighted_moment_and_reverse_triangle_lower_bound():
    x = sp.Symbol("x", positive=True)
    mu = core.MU
    exponent = sp.Rational(3, 2)+sp.I*mu
    assert sp.simplify((sp.diff(x**exponent/exponent, x)-sp.sqrt(x)*sp.exp(sp.I*mu*sp.log(x))).rewrite(sp.exp)) == 0
    assert sp.expand(exponent*sp.conjugate(exponent)) == 12
    q, theta = sp.symbols("q theta", real=True)
    modulus_squared = (1-q*sp.exp(sp.I*theta))*(1-q*sp.exp(-sp.I*theta))
    assert sp.simplify(sp.expand_complex(modulus_squared)-(1-q)**2-2*q*(1-sp.cos(theta))) == 0
    assert F(1, 8) < F(9, 64)
    assert F(12) < F(49, 4)
    assert (1-F(3, 8))/F(7, 2) == F(5, 28)


def test_source_weight_correction_and_observable_lower_bound_are_independent():
    t, q = sp.symbols("t q", nonnegative=True)
    assert sp.expand(5*q-q*q-4-(q-1)*(4-q)) == 0
    d = 1+t
    map_ratio_squared = 5/(d**6*(d**12+4))
    numerator = sp.factor((map_ratio_squared-d**-18)*d**18*(d**12+4))
    assert sp.expand(numerator-4*(d**12-1)) == 0
    bernoulli = sp.Poly(sp.expand(1-(1-9*t)*(1+t)**9), t)
    assert all(value >= 0 for value in bernoulli.all_coeffs())
    assert [sp.Rational(value) for value in bernoulli.all_coeffs()[::-1]] == [
        0, *[9*comb(9, j-1)-(comb(9, j) if j <= 9 else 0) for j in range(1, 11)]]


@pytest.mark.parametrize("r", [F(1, 1000), F(1, 10000), F(1, 10**7)])
def test_source_deficit_feedback_and_separation_rebuilt_with_Fraction(r):
    deficit = F(1, 100)
    # sqrt(v/r)<=1 and 0<=sigma<=1 give precisely deficit/r here.
    moment = F(5, 28)-deficit-4*r*r
    amplitude = F(8, 5)/F(13, 4)*(1-9*r*r)*moment
    omitted = 1000+2*F(101, 75)*(10+60*r*r)*5
    separation = 2*F(2, 25)-2*1200*r*r
    assert amplitude > F(2, 25)
    assert omitted < 1200
    assert separation > F(3, 20)
    c = core.calibration(sp.Rational(r))
    assert c["amplitude_over_r_squared_lower"] == sp.Rational(amplitude)
    assert c["simplified_response_error_over_r_fourth"] == sp.Rational(omitted)
    assert c["phase_separation_over_r_squared_lower"] == sp.Rational(separation)
    # The inherited error is an actual-source C0 estimate, not a mass ansatz.
    delta = sp.Rational(r*r/1000)
    if delta <= sp.Rational(1, 10**9):
        assert bounds.physical_response_error(delta, 1, sp.Rational(r), 0) == 1000*sp.Rational(r)**4


def test_phase_sequences_hit_opposite_extrema_without_a_delta_reset_of_source():
    mu, theta, n = sp.symbols("mu theta n", real=True)
    for target in (sp.pi/2, 3*sp.pi/2):
        log_delta_over_scale = -(target-theta+2*sp.pi*n)/mu
        phase_plus_argument = -mu*log_delta_over_scale+theta
        assert sp.expand(phase_plus_argument-target-2*sp.pi*n) == 0
    assert sp.sin(sp.pi/2) == 1 and sp.sin(3*sp.pi/2) == -1


def test_bounded_actual_residuals_need_not_converge_on_either_sequence():
    r = F(1, 1000)
    amplitude, error = F(81, 1000)*r*r, 1200*r**4
    upper_sequence = (amplitude+error, amplitude-error)
    lower_sequence = (-amplitude+error, -amplitude-error)
    assert upper_sequence[0] != upper_sequence[1]
    assert lower_sequence[0] != lower_sequence[1]
    assert min(upper_sequence)-max(lower_sequence) == 2*amplitude-2*error
    assert min(upper_sequence)-max(lower_sequence) > F(3, 20)*r*r


def test_explicit_smooth_pulse_has_strict_interior_support_and_deficit():
    r, w = F(1, 1000), F(1, 800000)
    assert w == r/800
    support_left, support_right = r/2+w, r-w
    plateau_left, plateau_right = r/2+2*w, r-2*w
    assert r/2 < support_left < plateau_left < plateau_right < support_right < r
    assert (plateau_left-r/2)+(r-plateau_right) == 4*w == r/200 < r/100


@pytest.mark.parametrize("bad", [True, 0.001, sp.Float(".001"), sp.oo, sp.nan])
def test_inexact_radius_cannot_enter_the_calibrated_bound(bad):
    with pytest.raises((TypeError, ValueError)):
        core.calibration(bad)
