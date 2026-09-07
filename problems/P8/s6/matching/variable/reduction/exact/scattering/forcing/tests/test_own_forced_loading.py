"""Exact fixed pulse, literal canonical source, and positive loading margins."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_own_forced import loading as l


def test_primary_exact_identities():
    assert l.identities()
    assert all(value == 0 for value in l.identities().values())


def test_all_primary_margins():
    assert all(l.checks().values())


def test_independent_fraction_source_bounds():
    f = Fraction
    length, delta = f(1, 100), f(1, 10**6)
    pole_floor = 16 / (delta + 8 * length**2)
    potential_floor = pole_floor - 44
    mass_floor = potential_floor - 22
    kinetic_floor = f(19, 10)**3 / f(21, 10)
    kinetic_ceiling = f(21, 10)**3 / f(19, 10)
    assert mass_floor == f(15947134, 801) > 1 / length**2
    assert potential_floor > 0
    assert kinetic_floor == f(6859, 2100) > 1
    assert kinetic_ceiling == f(9261, 1900) < 5
    actual = l.source_bounds()
    assert actual["mass_lower"] == mass_floor
    assert actual["kinetic_lower"] == kinetic_floor
    assert actual["kinetic_upper"] == kinetic_ceiling


def test_source_retains_pump_and_sqrt_kinetic():
    kinetic, mass, pump, q = sp.symbols("k mass pump q", positive=True)
    y, yy = sp.symbols("Y Y_uu")
    full = yy + (mass - pump) * y - sp.sqrt(kinetic) * mass * q
    omitted = yy + (mass - pump) * y - sp.sqrt(kinetic) * (mass - pump) * q
    assert sp.expand(full - omitted) == -sp.sqrt(kinetic) * pump * q


def test_independent_green_first_zero_constants():
    f = Fraction
    length = f(1, 100)
    potential = 8 / length**2 + 44
    width = length / 4
    derivative = 1 - potential * width**2 / 2
    ratio = 1 - potential * width**2 / 6
    assert derivative == f(59989, 80000) > f(2, 3)
    assert ratio == f(219989, 240000) > f(9, 10)
    assert l.volterra_bounds(potential, width)["G_over_lag_lower"] == ratio


def test_uniform_sharper_loading_is_kept_before_the_limit():
    data = l.calibration()
    assert data["Y_load_sharper_lower_over_eta"] == Fraction(219989, 30720000)
    assert data["Y_load_sharper_lower_over_eta"] > data["Y_load_lower_over_eta"]
    assert data["Y_load_lower_over_eta"] == Fraction(9, 1280)
    assert data["Y_u_load_lower_over_eta"] == Fraction(25, 6)


def test_pulse_is_fixed_and_has_strict_support_margins():
    p = l.pulse()
    assert p["support_left"] + 3 * l.L / 4 == l.L / 128
    assert -5 * l.L / 8 - p["support_right"] == l.L / 128
    assert p["source_to_load_lag_min"] == 17 * l.L / 128
    assert p["source_to_load_lag_max"] == 31 * l.L / 128
    assert p["plateau_width"] == 5 * l.L / 64 > l.L / 16
    assert p["q"].free_symbols == {l.U, l.ETA}


@pytest.mark.parametrize("point", [-l.L, -3 * l.L / 4, -5 * l.L / 8, -l.A, 0, l.A])
def test_pulse_zero_at_and_beyond_support(point):
    assert sp.simplify(l.pulse()["q"].subs(l.U, sp.Rational(point))) == 0


def test_plateau_height_is_exact():
    p = l.pulse()
    midpoint = (p["plateau_left"] + p["plateau_right"]) / 2
    assert sp.simplify(p["q"].subs(l.U, sp.Rational(midpoint))) == l.ETA


@pytest.mark.parametrize("order", range(7))
def test_flat_endpoint_derivative_recurrence(order):
    x = sp.Symbol("x", positive=True)
    z = sp.Symbol("z", real=True)
    actual = sp.diff(sp.exp(-1 / x), x, order) / sp.exp(-1 / x)
    expected = l.flat_derivative_polynomial(order).subs(z, 1 / x)
    assert sp.simplify(actual - expected) == 0


@pytest.mark.parametrize("bad", [True, False, 1.0, sp.Float(1), -1])
def test_flat_derivative_order_guard_after_warming(bad):
    l.flat_derivative_polynomial(1)
    with pytest.raises((TypeError, ValueError)):
        l.flat_derivative_polynomial(bad)


def test_punctured_limit_has_a_denominator_margin():
    f = Fraction
    a = f(1, 200)
    lower = 2 - 2 / (1 + a**2)**4
    assert lower > 0
    assert l.source_bounds()["punctured_D_lower"] == lower
    # The center itself has no such delta=0 literal-coefficient margin.
    assert 2 - 2 / (1 + f(0))**4 == 0


def test_too_long_first_zero_interval_is_only_an_omission_control():
    failed = l.volterra_bounds(80044, l.L / 2)
    assert failed["positive_derivative_bootstrap"] is False
    assert failed["derivative_lower"] < 0


@pytest.mark.parametrize("bad", [True, 1.0, sp.oo, sp.nan])
def test_exact_bound_input_guards(bad):
    with pytest.raises((TypeError, ValueError)):
        l.volterra_bounds(bad, l.L / 4)


def test_negative_bound_domain_is_rejected():
    with pytest.raises(ValueError):
        l.volterra_bounds(-1, l.L / 4)
    with pytest.raises(ValueError):
        l.volterra_bounds(1, -1)
