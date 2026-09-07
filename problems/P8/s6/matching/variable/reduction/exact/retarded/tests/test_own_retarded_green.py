from fractions import Fraction

import pytest
import sympy as sp
from p8_own_retarded import green


def rational(value):
    return sp.Rational(value.numerator, value.denominator)


def test_all_exact_margins():
    assert len(green.identities()) == 14
    assert set(green.identities().values()) == {0}
    assert len(green.checks()) == 22
    assert all(green.checks().values())


def test_independent_positive_flux_integrals():
    h, t = sp.symbols("h t", nonnegative=True)
    kmin, kmax, smax = map(rational, (green.K_MIN, green.K_MAX, green.S_MAX))
    flux = 1-sp.integrate(smax*t/kmin, (t, 0, h))
    lower = sp.integrate((1-smax*t**2/(2*kmin))/kmax, (t, 0, h))
    endpoint = {h: sp.Rational(1, 4)}
    assert flux.subs(endpoint) == sp.Rational(283, 608)
    assert (flux/kmax).subs(endpoint) == sp.Rational(1415, 12768)
    assert sp.cancel(lower/h).subs(endpoint) == sp.Rational(7495, 38304)
    assert sp.diff(flux, h) == -smax*h/kmin
    assert sp.diff(sp.cancel(lower/h), h) == -smax*h/(3*kmin*kmax)


def test_direct_volterra_route_is_weaker_but_sufficient():
    h = sp.Rational(1, 4)
    kmin, kmax, smax = map(rational, (green.K_MIN, green.K_MAX, green.S_MAX))
    direct = 1/kmax-smax*h**2/(6*kmin**2)
    sharp = (1-smax*h**2/(6*kmin))/kmax
    assert direct == sp.Rational(46385, 242592)
    assert sp.Rational(3, 16) < direct < sharp


def test_full_physical_time_and_source_normalization():
    x = green.X
    M, tau, delta = sp.symbols("M tau delta", positive=True)
    k, hidden, prescribed = (sp.Function(name)(x) for name in ("k", "Q", "q"))
    numerator, b, Dbar = sp.symbols("profile_Q b D_over_delta", positive=True)
    beta1 = M**2*numerator/(tau**2*delta*Dbar)
    source_coefficient = 2*numerator*b/Dbar
    physical_D = lambda value: sp.diff(value, x)/(tau*sp.sqrt(delta))
    physical = physical_D(M**2*k*physical_D(hidden))+2*beta1*b*(hidden-prescribed)
    normalized = sp.diff(k*sp.diff(hidden, x), x)+source_coefficient*(hidden-prescribed)
    assert sp.simplify(physical*tau**2*delta/M**2-normalized) == 0
    assert source_coefficient != numerator*b/Dbar


def test_literal_sturm_to_volterra_identity():
    x, a, t, z = sp.symbols("x a t z", real=True)
    k, s, R = sp.Function("k"), sp.Function("s"), sp.Function("R")
    flux = 1-sp.Integral(s(z)*R(z), (z, a, x))
    integral_R = sp.Integral((1-sp.Integral(s(z)*R(z), (z, a, t)))/k(t), (t, a, x))
    assert sp.simplify(k(x)*sp.diff(integral_R, x)-flux) == 0
    assert sp.diff(flux, x)+s(x)*R(x) == 0
    assert integral_R.subs(x, a).doit() == 0
    assert flux.subs(x, a).doit() == 1


def test_constant_coefficient_green_and_wrong_jump_control():
    x, a = sp.symbols("x a", real=True)
    R = sp.sin(4*(x-a))/16
    assert sp.simplify(sp.diff(R, x, 2)*4+64*R) == 0
    assert R.subs(x, a) == 0
    assert (4*sp.diff(R, x)).subs(x, a) == 1
    wrong = sp.sin(4*(x-a))/4
    assert (4*sp.diff(wrong, x)).subs(x, a) == 4
    # Independent Taylor inequalities on 0<=4h<=1 give stronger toy bounds.
    assert sp.Rational(1, 4)*(1-sp.Rational(1, 2)) > sp.Rational(1, 10)
    assert sp.Rational(1, 4)-sp.Rational(2, 3)/16 > sp.Rational(3, 16)


def test_time_dependent_kinetic_green_normalization():
    x, a = sp.symbols("x a", real=True)
    alpha = sp.Rational(1, 10)
    k = 4*sp.exp(alpha*x)
    s = k*(16+alpha**2/4)
    R = sp.exp(-alpha*(x+a)/2)*sp.sin(4*(x-a))/16
    assert sp.simplify(sp.diff(k*sp.diff(R, x), x)+s*R) == 0
    assert sp.simplify((k*sp.diff(R, x)).subs(x, a)) == 1
    # Dropping k' changes the same literal differential equation.
    assert sp.simplify(k*sp.diff(R, x, 2)+s*R) != 0


def test_pulse_lower_bounds_and_metric_amplitude_upper_bounds():
    c = green.response_constants()
    minimum_distance = Fraction(1, 8)
    area = Fraction(1, 16)
    assert Fraction(3, 16)*minimum_distance*42*area == Fraction(63, 1024)
    assert Fraction(1, 10)*42*area == Fraction(21, 80)
    assert c["Q_final_sharper_lower_over_eta"] > c["Q_final_lower_over_eta"]
    assert c["Q_x_final_sharper_lower_over_eta"] > c["Q_x_final_lower_over_eta"]
    assert c["Q_whole_window_upper_over_eta"] == Fraction(975, 2432) < 1
    assert c["Q_x_whole_window_upper_over_eta"] == Fraction(325, 152)


def test_pulse_class_has_exact_support_and_sufficient_plateau():
    p = green.pulse()
    assert green.LEFT < p["support_left"] < p["plateau_left"]
    assert p["plateau_right"] < p["support_right"] < green.PULSE_RIGHT
    assert p["plateau_right"]-p["plateau_left"] == Fraction(5, 64)
    assert p["area_margin_over_eta"] == Fraction(1, 64)
    for point in (Fraction(0), Fraction(-1, 128), green.LEFT, p["support_left"], p["support_right"]):
        assert p["profile"].subs(green.X, rational(point)) == 0
    midpoint = (p["plateau_left"]+p["plateau_right"])/2
    assert p["profile"].subs(green.X, rational(midpoint)) == green.ETA


@pytest.mark.parametrize("order", range(7))
def test_flat_step_derivative_recurrence(order):
    data = green.flat_derivative_polynomials(order)
    t = sp.Symbol("t", positive=True)
    polynomial = data["polynomials"][order]
    derivative = sp.diff(sp.exp(-1/t), t, order)
    expected = sp.exp(-1/t)*polynomial.as_expr().subs(data["z"], 1/t)
    assert sp.simplify(derivative-expected) == 0
    assert polynomial.degree() == 2*order
    assert sp.limit(sp.exp(-1/t)*polynomial.as_expr().subs(data["z"], 1/t), t, 0, dir="+") == 0


def test_same_final_germ_and_distinct_retarded_history():
    p = green.pulse()
    assert p["support_right"] < green.PULSE_RIGHT < 0
    # The pulse is identically zero on the OPEN interval (support_right,+infty).
    # Thus its entire final germ equals that of zero input, not only finite jets.
    assert p["profile"].subs(green.X, 0) == 0
    response = green.response_constants()
    assert response["Q_final_lower_over_eta"] > 0
    assert response["Q_x_final_lower_over_eta"] > 0
    # Any common output for these two germs incurs at least half the gap.
    assert response["Q_final_lower_over_eta"]/2 == Fraction(63, 2048)


def test_physical_clock_derivative_is_not_fixed_frequency():
    response = green.response_constants()
    assert response["physical_Q_time_derivative_lower"] == 21*green.ETA/(80*green.TAU*sp.sqrt(green.DELTA))
    assert response["physical_window_duration"] == green.TAU*sp.sqrt(green.DELTA)/4
    assert sp.limit(response["physical_window_duration"], green.DELTA, 0, dir="+") == 0
    assert sp.limit(response["physical_Q_time_derivative_lower"], green.DELTA, 0, dir="+") is sp.oo


def test_longer_window_is_not_silently_certified():
    failed = green.volterra_bounds(width=Fraction(1, 2))
    assert failed["flux_lower"] < 0
    assert not failed["positive_flux_bootstrap"]
    assert green.volterra_bounds()["positive_flux_bootstrap"]


@pytest.mark.parametrize("bad", [True, 0.25, "1/4", float("inf"), float("nan")])
def test_exact_bound_inputs(bad):
    with pytest.raises(TypeError):
        green.volterra_bounds(width=bad)


@pytest.mark.parametrize("parameters", [{"k_min": 0}, {"k_min": 5, "k_max": 4}, {"s_max": -1}, {"width": -1}])
def test_invalid_bound_domains(parameters):
    with pytest.raises(ValueError):
        green.volterra_bounds(**parameters)


def test_derivative_order_validation_precedes_computation():
    assert len(green.flat_derivative_polynomials(1)["polynomials"]) == 2
    for bad in (True, 1.0, sp.Float(1), "1"):
        with pytest.raises(TypeError):
            green.flat_derivative_polynomials(bad)
    with pytest.raises(ValueError):
        green.flat_derivative_polynomials(-1)
