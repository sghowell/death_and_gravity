"""Independent phase, boundary, radial and proper-normalization audits."""

import sympy as sp
from p8a_existence import mode_lipschitz as modes


def test_actual_two_mode_quadratic_kernel_not_only_constant_potential():
    k, t, s, r = sp.symbols("k t s r", positive=True)
    # Ordered domain r<s<t: literal w2*w0bar + c.c. and the two halves
    # of w1*w1bar. Neither individual contribution is the full response.
    ordered_w2 = 2*sp.sin(k*(t-s))*sp.sin(k*(s-r))*sp.cos(k*(t-r))/k**2
    ordered_w1_squared = 2*sp.sin(k*(t-s))*sp.sin(k*(t-r))*sp.cos(k*(s-r))/k**2
    full = (sp.cos(2*k*(s-r))-sp.cos(2*k*(t-r)))/k**2
    assert sp.trigsimp(sp.expand_trig(ordered_w2+ordered_w1_squared-full)) == 0
    assert sp.trigsimp((ordered_w2-full).subs({t: 3, s: 2, r: 1, k: sp.pi/4})) != 0


def test_abel_quadratic_kernel_from_frequency_parameter_differentiation():
    a, b, epsilon = sp.symbols("a b epsilon", positive=True)
    radial = sp.log((epsilon**2+4*b**2)/(epsilon**2+4*a**2))/2
    # Differentiate the convergent Abel radial integral wrt each endpoint;
    # integral exp(-epsilon*k) sin(2*b*k) dk = 2b/(epsilon^2+4b^2).
    assert sp.diff(radial, b) == 4*b/(epsilon**2+4*b**2)
    assert sp.diff(radial, a) == -4*a/(epsilon**2+4*a**2)
    assert sp.simplify(radial.subs(b, a)) == 0
    assert sp.simplify(sp.limit(radial, epsilon, 0, dir="+")-sp.log(b/a)) == 0


def test_abel_error_bounds_have_correct_half_frequency_factors():
    x, epsilon = sp.symbols("x epsilon", positive=True)
    derivative_error_mass = sp.integrate((epsilon/2)**2/(x**2+(epsilon/2)**2), (x, 0, sp.oo))
    assert sp.simplify(derivative_error_mass-sp.pi*epsilon/4) == 0
    log_error_mass = sp.integrate(sp.log(1+(epsilon/(2*x))**2), (x, 0, sp.oo))
    assert sp.simplify(log_error_mass-sp.pi*epsilon/2) == 0


def test_derivative_transfer_retains_the_nonflat_past_endpoint():
    t, s, k = sp.symbols("t s k", positive=True)
    u0, u1 = sp.symbols("u0 u1", real=True)
    potential = u0+u1*s
    green = sp.sin(k*(t-s))/k
    wave = sp.exp(-sp.I*k*s)
    # The integral kernel of w1'+ik*w1-Dw1[U](U') is a total s derivative.
    integrand = -(sp.cos(k*(t-s))+sp.I*sp.sin(k*(t-s)))*potential*wave+green*u1*wave
    boundary_primitive = green*potential*wave
    assert sp.simplify(integrand-sp.diff(boundary_primitive, s)) == 0
    assert sp.simplify(boundary_primitive.subs(s, t)-boundary_primitive.subs(s, 0)+sp.sin(k*t)*u0/k) == 0
    assert -sp.sin(k*t)*u0/k != 0


def test_shared_history_quadratic_bound_is_derived_not_assumed():
    t, length, v = sp.symbols("t L v", positive=True)
    outer_mark = sp.integrate((length-v)*sp.log(t/v), (v, 0, length))
    inner_mark = length**2/2
    assert sp.simplify(outer_mark+inner_mark-length**2*(sp.Rational(5, 4)+sp.log(t/length)/2)) == 0
    cap = length**2*(sp.Rational(5, 4)+sp.log(t/length)/2)
    assert sp.simplify(sp.diff(cap, length)-length*(2+sp.log(t/length))) == 0


def test_quadratic_radial_measure_and_oscillator_normalization():
    k, t = sp.symbols("k t", positive=True)
    coefficient = sp.sin(k*t)**2/k**4-t*sp.sin(2*k*t)/(2*k**3)
    primitive = -sp.sin(k*t)**2/(2*k**2)
    assert sp.trigsimp(sp.diff(primitive, k)-k*coefficient) == 0
    assert sp.limit(primitive, k, sp.oo)-sp.limit(primitive, k, 0, dir="+") == t**2/2


def test_actual_Wick_response_clock_and_hbar_normalization():
    scale, hbar, a, derivative = sp.symbols("eta_star hbar a R_x", positive=True)
    # x=(eta-eta_start)/eta_star and physical U=U_dimensionless/eta_star^2.
    # k dk contributes eta_star^-2 to Rmode; one eta derivative adds -1.
    quantum_response = hbar/(4*sp.pi**2*a**2*scale**2)
    clock_response = quantum_response*derivative/scale
    assert sp.simplify(clock_response-hbar*derivative/(4*sp.pi**2*a**2*scale**3)) == 0
    assert modes.PREPARED_U1 == sp.Rational(61013499, 8192)
