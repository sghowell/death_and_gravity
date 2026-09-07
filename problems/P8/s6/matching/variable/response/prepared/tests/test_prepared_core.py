import pytest
import sympy as sp
from p8_variable_prepared import analytic, exact, independent


def test_full_parent_analytic_residuals_are_exact():
    residuals = analytic.checks()
    assert len(residuals) == 23
    assert all(sp.simplify(value) == 0 for value in residuals.values())
    assert not any(value.has(sp.Float) for value in residuals.values())


@pytest.mark.parametrize("j,n", [(0, 0), (0, 1), (1, 2), (0, 7), (3, 12)])
def test_inverse_matches_closed_fraction_product(j, n):
    d = analytic.derive()
    polynomial = sp.Poly(analytic.inverse_monomial(j, n), d["delta"], d["u"])
    actual = {powers: sp.Rational(value) for powers, value in polynomial.terms()}
    assert actual == {powers: sp.Rational(value) for powers, value in independent.inverse_monomial(j, n).items()}


def test_omitting_delta_shift_is_not_an_inverse():
    d = analytic.derive()
    wrong = d["u"]**2/analytic.diagonal(2)
    assert sp.factor(analytic.leading_operator(wrong)-d["u"]**2) == d["delta"]/88


def test_moving_weight_omission_changes_even_quartic():
    d = analytic.derive()
    momentum, u, delta = d["K"], d["u"], d["delta"]
    forcing_without_f = -2*momentum*delta/5-(16*momentum/5-sp.Rational(864, 25))*u**2
    full = analytic.leading_operator(analytic.low_jets()["q_even"])
    assert sp.expand(full-forcing_without_f+48*(5*momentum+18)*u**2/25) == 0


def test_K_zero_is_an_exact_common_mode_not_a_division_assumption():
    c = analytic.checks()
    assert c["K_zero_common_light"] == c["K_zero_common_heavy"] == 0
    assert analytic.low_jets()["q_even"].subs(analytic.derive()["K"], 0) == 0


@pytest.mark.parametrize("bad", [True, 0.1, sp.Float("0.1"), sp.oo, sp.nan, sp.I, sp.Symbol("x")])
def test_exact_real_domain_rejects_unproved_values(bad):
    with pytest.raises((TypeError, ValueError)):
        exact.number(bad)


@pytest.mark.parametrize("bad", [0, -1, sp.Rational(1, 10**8), True])
def test_actual_delta_domain_is_not_the_complex_proof_disc(bad):
    with pytest.raises((TypeError, ValueError)):
        exact.parameters(bad)


def test_K_zero_new_analytic_domain_but_not_old_transfer_domain():
    delta = sp.Rational(1, 10**10)
    assert exact.parameters(delta, 0) == (delta, 0)
    with pytest.raises(ValueError):
        exact.parameters(delta, 0, transfer=True)
    with pytest.raises(ValueError):
        exact.parameters(delta, sp.Rational(401, 100))
