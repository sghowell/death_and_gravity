"""Exact index endpoints, genuine history condition and relative-pole moments."""

import pytest
import sympy as sp
from p8a_maxwell_focusing import history


def test_every_history_index_and_sampler_identity_is_zero():
    assert all(sp.simplify(value) == 0 for value in history.identities().values())


def test_cubic_join_has_the_two_required_H2_traces():
    x = sp.Symbol("x")
    p = history.cubic(x)
    assert [p.subs(x, 0), sp.diff(p, x).subs(x, 0)] == [0, 0]
    assert [p.subs(x, 1), sp.diff(p, x).subs(x, 1)] == [1, 0]


def test_instantaneous_K_does_not_supply_a_quantitative_history_duration():
    # H(0) can be fixed while the earlier squared history changes.
    x = sp.Symbol("x", real=True)
    u = history.cubic(x)
    h1, h2 = -sp.Rational(3, 2), -sp.Rational(3, 2)+2*(1-x)
    gains = [sp.integrate(3*((sp.diff(u, x)-h*u)**2-sp.diff(u, x)**2), (x, 0, 1))
             for h in (h1, h2)]
    assert h1 == h2.subs(x, 1)
    assert gains[0] != gains[1]


def test_history_bound_requires_positive_length_and_contraction():
    with pytest.raises(ValueError):
        history.lower(1, 0)
    with pytest.raises(ValueError):
        history.lower(0, 1)


def test_linear_endpoint_cannot_replace_the_cubic_in_relative_pole_cost():
    x = sp.Symbol("x", positive=True)
    assert sp.integrate((x/x**2)**2, (x, 0, 1)) == sp.oo
    assert history.moments()["pole_zero"] == sp.Rational(13, 3)
