"""Exact nonlinear source identities and continuous Taylor controls."""
import pytest
import sympy as sp
from p8_aligned_elimination import bounds, source


def test_exact_ODE_and_full_source_remainders():
    assert all(value == 0 for value in source.checks().values())
    assert all(value is True for value in bounds.proof_checks().values())


@pytest.mark.parametrize("r", (-sp.Rational(1, 25), -sp.Rational(1, 100), 0,
                               sp.Rational(1, 100), sp.Rational(1, 25)))
@pytest.mark.parametrize("h,H", ((1, 0), (2, -2), (3, 2)))
def test_independent_full_source_values(r, h, H):
    r, h, H = map(sp.sympify, (r, h, H))
    data = source.expansion()
    f1 = -9*sp.Rational(H, 8*h**2)
    N = 1+r
    y = 1-N**-2
    # A generic exact cubic remainder checks the whole identity, not the
    # original ODE's bound by numerical samples.
    Q = f1*y**2/2+sp.Rational(2, 3)*y**3
    k = sp.Rational(1, 100)
    original = (N**-2-1)*(k+3*H*(1-1/N))/h+sp.Rational(3, 2)*Q/N
    substitution = {source.r: r, source.h: h, source.H: H, source.k: k,
                    source.f1: f1, source.RQ: sp.Rational(2, 3)*y**3}
    assert sp.factor(original-(data["quadratic"]+data["normal_remainder"]).subs(substitution)) == 0
    assert sp.factor(N*original-(data["quadratic"]+data["coordinate_remainder"]).subs(substitution)) == 0


def test_derivative_control_keeps_Q_remainder_chain_rule():
    data = source.expansion()
    wrong = sp.diff(data["coordinate_remainder"], source.r)
    assert sp.factor(data["coordinate_r_remainder"]-wrong) == 3*source.RQx/(1+source.r)**3


def test_first_nonlinear_remainder_is_not_identically_zero():
    data = source.expansion()
    difference = data["normal_remainder"].subs({source.H: 0, source.f1: 0, source.RQ: 0})
    assert sp.factor(difference) != 0
    assert sp.limit(difference/source.r**2, source.r, 0) == 3*source.k/source.h


def test_normalized_source_and_density_examples():
    result = bounds.values(sp.Rational(1, 1000), sp.Rational(1, 10**9))
    assert result["normal_source_remainder"] == sp.Rational(13, 125000000)
    assert result["physical_electric_source_remainder"] == sp.Rational(3, 10**7)
    assert result["local_first_derivative_Maxwell_density_remainder"] == sp.Rational(573, 4*10**22)
    assert result["nonlocal_inverse_or_nonlinear_solution_bound"] is False


def test_nonunit_normalization_is_not_lost():
    result = bounds.physical_values(sp.Rational(1, 1000), sp.Rational(1, 10000), 3, 2)
    assert result["normalized_zeta"] == sp.Rational(1, 120000)
    assert result["normal_source_remainder"] == sp.Rational(13, 250000000)
    assert result["electric_source_remainder"] == sp.Rational(3, 4*10**7)
    assert result["local_Maxwell_density_remainder"] == sp.Rational(573, 64*10**17)


def test_preparation_and_nonzero_third_source_survive_the_full_source_map():
    data = source.preparation()
    assert data["prepared_full_source_jet_through_three"] == 0
    assert data["third_source_jet"] != 0


def test_full_lapse_map_is_not_a_finite_Fourier_polynomial():
    # Even its coefficient of k is rational in r=epsilon*cos(theta).
    data = source.expansion()
    coefficient = sp.factor(sp.diff(data["coordinate"], source.k))
    assert sp.denom(coefficient).has(1+source.r)
    assert sp.factor(sp.diff(coefficient, source.r, 7)) != 0


@pytest.mark.parametrize("epsilon", (0, -1, sp.Rational(1, 20)))
def test_outside_lapse_domain_is_rejected(epsilon):
    with pytest.raises(ValueError):
        bounds.values(epsilon, 1)
