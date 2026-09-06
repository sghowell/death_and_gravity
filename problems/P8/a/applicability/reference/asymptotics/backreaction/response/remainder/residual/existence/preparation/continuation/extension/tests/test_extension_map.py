"""Retained actual source terms and original-data exclusion controls."""

import sympy as sp
from p8a_extension import map as equation


def test_full_forced_source_and_map_identities():
    assert all(value == 0 for value in equation.identities().values())


def test_source_difference_is_not_silently_omitted_after_resummation():
    q = sp.symbols("q0:7")
    result = equation.full_forced_remainder(*q)
    assert sp.diff(result, q[5]) == 1
    assert result.subs({q[0]: 0, q[1]: 0, q[2]: 0, q[3]: 0, q[4]: 0, q[5]: 1}) == 1


def test_both_source_endpoints_and_hubble_dependence_remain():
    c, h, initial, history = sp.symbols("c h initial history", positive=True)
    source = equation.source_primitive(c, h, initial, history)
    assert sp.diff(source, initial) == -8
    assert sp.diff(source, history) == -8
    assert sp.diff(source, h) == -8*c/h**2
    assert source.subs({initial: c/h, history: 0}) == 0


def test_switching_density_with_radiation_pressure_is_a_different_source():
    t = sp.Symbol("t", real=True)
    a, c = sp.Function("a", positive=True)(t), sp.Function("c")(t)
    h = sp.diff(a, t)/a
    rho, pressure = c/a**4, c/(3*a**4)
    assert sp.simplify(sp.diff(rho, t)+3*h*(rho+pressure)) == sp.diff(c, t)/a**4
