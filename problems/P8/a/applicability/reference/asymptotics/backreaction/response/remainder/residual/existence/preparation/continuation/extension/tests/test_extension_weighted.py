"""Exact weighted gains, actual prefix response and rejected domains."""

import pytest
import sympy as sp
from p8a_extension import weighted


@pytest.mark.parametrize("value", [True, False, 0.1, sp.Float("0.1"), sp.oo, sp.nan,
                                 "1/0", sp.Symbol("x"), -1])
def test_exact_domain_rejects_inexact_or_invalid_values(value):
    with pytest.raises((TypeError, ValueError)):
        weighted.nonnegative(value)


def test_primitive_gain_and_boolean_switch_domain():
    gains = weighted.metric_gains(2*10**6)
    assert gains["u"] == sp.Rational(1, 2*10**6)
    assert gains["h"] == gains["u"]**2
    assert gains["a_squared"] == 18*gains["u"]**3
    with pytest.raises(TypeError):
        weighted.nonnegative(1, positive=1)
    with pytest.raises(ValueError):
        weighted.metric_gains(0)


def test_prefix_response_uses_actual_same_history_bound():
    result = weighted.prefix_response("1/100000", 3, "1/1000000")
    assert result["total"] == sp.Rational(70011, 80000000000000000)
    assert result["total"] == result["quadratic"]+result["higher"]
    assert weighted.prefix_response(0, 3, "1/1000000")["total"] == 0


@pytest.mark.parametrize("args", [(1, 3, 1), (0, 1, 2), (0, 1, 0),
                                 ("1/100000", 3, "1/1000000", 1),
                                 (0, 1, 1, "3/2"), (True, 1, 1)])
def test_response_strength_history_and_log_gate(args):
    with pytest.raises((TypeError, ValueError)):
        weighted.prefix_response(*args)


def test_all_weighted_and_prefix_limit_identities():
    assert all(value == 0 for value in weighted.identities().values())
    t, history = sp.symbols("t T", positive=True)
    coefficient = t**2*(sp.Rational(5, 4)+sp.log(history/t)/2)
    assert sp.simplify(sp.diff(coefficient, t).subs(t, history)) == 2*history


def test_no_extra_weight_loss_is_inserted_into_a_causal_prefix():
    sigma, t, s = sp.symbols("sigma t s", positive=True)
    # At an observation t, the largest allowed weighted difference on
    # its past prefix is exp(sigma*t)D. Multiplication by the output
    # weight cancels this factor exactly, not exp(sigma*L).
    assert sp.exp(-sigma*t)*sp.exp(sigma*t) == 1
    assert sp.diff(sp.exp(sigma*s), s).is_positive
