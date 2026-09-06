"""Sheet, sign, normalized inverse and forced-response controls."""

import pytest
import sympy as sp
from p8a_continuation import resolvent


def test_all_formal_pole_cut_and_scaling_identities():
    assert all(value == 0 for value in resolvent.identities().values())


def test_principal_sheet_pole_list_keeps_the_damped_pair():
    beta = sp.Symbol("beta", real=True)
    c = sp.Symbol("c", positive=True)
    poles = resolvent.pole_expressions(beta, c)
    assert set(poles) == {0, 1, -1}
    assert poles[0] == sp.exp(sp.LambertW(2*c*sp.exp(2*beta))/2-beta)
    assert poles[1] != poles[0]
    assert poles[-1] != poles[0]


def test_positive_axis_monotonicity_and_unique_root_signs():
    s, c = sp.symbols("s c", positive=True)
    beta = sp.Symbol("beta", real=True)
    symbol = sp.log(s)+beta-c/s**2
    assert sp.diff(symbol, s) == 1/s+2*c/s**3
    assert sp.diff(symbol, s).is_positive
    assert sp.simplify(symbol.subs(s, sp.exp(-beta))) == -c*sp.exp(2*beta)
    assert sp.limit(symbol, s, 0, dir="+") == -sp.oo
    assert sp.limit(symbol, s, sp.oo) == sp.oo


def test_wrong_Einstein_sign_has_a_different_minimum():
    s, c = sp.symbols("s c", positive=True)
    beta = sp.Symbol("beta", real=True)
    wrong = sp.log(s)+beta+c/s**2
    assert sp.simplify(sp.diff(wrong, s).subs(s, sp.sqrt(2*c))) == 0
    assert sp.simplify(wrong.subs(s, sp.sqrt(2*c))-(sp.log(2*c)+2*beta+1)/2) == 0
    correct = 1/resolvent.transfer(s, beta, c)
    assert sp.simplify(wrong/2-correct) == c/s**2


def test_cut_uses_full_Einstein_term_even_on_negative_axis():
    r, c = sp.symbols("r c", positive=True)
    beta = sp.Symbol("beta", real=True)
    correct = resolvent.cut_density(r, beta, c)
    omitted = resolvent.cut_density(r, beta, 0)
    assert correct != omitted
    assert sp.simplify(correct.subs(r, 1)) == 2/((beta-c)**2+sp.pi**2)
    assert correct.is_positive


def test_pole_removal_fails_the_inverse_identity():
    s, p, residue = sp.symbols("s p residue", positive=True)
    beta, c = sp.symbols("beta c", real=True)
    exact = resolvent.transfer(s, beta, c)
    block = (sp.log(s)+beta-c/s**2)/2
    assert sp.simplify(block*exact-1) == 0
    assert sp.simplify(block*(exact-residue/(s-p))-1+residue*block/(s-p)) == 0
    assert sp.simplify(-residue*block/(s-p)) != 0


def test_causal_pole_response_includes_forcing_history_not_free_amplitude():
    p, t, length, residue = sp.symbols("p t length residue", positive=True)
    s = sp.Symbol("s", nonnegative=True)
    while_on = sp.integrate(residue*sp.exp(p*(t-s)), (s, 0, t))
    after = sp.integrate(residue*sp.exp(p*(t-s)), (s, 0, length))
    assert sp.simplify(while_on-residue*(sp.exp(p*t)-1)/p) == 0
    assert sp.simplify(after-residue*sp.exp(p*t)*(1-sp.exp(-p*length))/p) == 0
    assert while_on.subs(t, 0) == 0
    assert sp.exp(p*t).subs(t, 0) != 0


@pytest.mark.parametrize("value", [True, 0.1, sp.Float("0.1"), sp.oo, sp.nan, "1/0", -1, 0])
def test_named_parameters_reject_inexact_or_nonpositive_scales(value):
    with pytest.raises((TypeError, ValueError)):
        resolvent.named_parameters(scale=value)
    with pytest.raises((TypeError, ValueError)):
        resolvent.named_parameters(delta=value)
