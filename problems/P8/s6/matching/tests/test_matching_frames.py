import pytest
import sympy as sp
from p8_matching import frames


def test_independent_tensor_pullback_and_inverse():
    result = frames.identities()
    assert all(sp.simplify(value) == 0 for value in result["residuals"].values())
    assert result["all_time_threshold_limit"] == 1


@pytest.mark.parametrize("C,W", [(1, 1), (4, 9), (sp.Rational(1, 4), sp.Rational(9, 4))])
def test_tensor_roundtrip_and_conformal_limit(C, W):
    g, f = frames.tensor_map(C, W)
    assert sp.simplify(sp.sqrt(g*f)-C) == 0
    assert sp.simplify(f**sp.Rational(3, 2)/sp.sqrt(g)-W) == 0
    assert frames.tensor_map(C, C) == (sp.sympify(C), sp.sympify(C))


def test_three_slice_obstruction_and_nonclaim_at_equality():
    assert frames.three_slice_bound(1, 0)["squared_scale_margin"] == "3"
    assert frames.three_slice_bound(1, sp.Rational(1, 10))["status"] == "EXCLUDED_BY_EINSTEIN_CONCAVITY"
    assert frames.three_slice_bound(1, sp.Rational(3, 5))["status"] == "BOUND_INCONCLUSIVE"
    assert frames.three_slice_bound(2, 0)["necessary_error_lower_bound"] == "12/13"
    assert frames.three_slice_bound(1, sp.Rational(9, 10))["status"] == "BOUND_INCONCLUSIVE"


@pytest.mark.parametrize("L,epsilon", [(0, 0), (-1, 0), (1, -1), (1, 1), (1, 0.1)])
def test_bad_slice_domains_rejected(L, epsilon):
    with pytest.raises(ValueError):
        frames.three_slice_bound(L, epsilon)


def test_nonminimal_escape_fails_two_tail_growth_not_concavity():
    u = sp.Symbol("u", real=True)
    q = sp.exp(-3*u**2)
    aE = sp.sqrt(q)*(1+u**2)
    assert sp.limit(aE, u, sp.oo) == sp.limit(aE, u, -sp.oo) == 0
    assert frames.tensor_map(q, q) == (q, q)
