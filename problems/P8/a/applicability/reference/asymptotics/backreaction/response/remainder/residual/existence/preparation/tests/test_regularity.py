import itertools
from math import factorial

import pytest
import sympy as sp
from p8a_preparation import regularity


def test_highest_jet_contraction_is_strict_on_actual_ball():
    data = regularity.highest_jet_bound(
        sp.Rational(1, 10**5), 3, sp.Rational(17, 30), sp.Rational(21, 100))
    assert data["strength"] == sp.Rational(27, 10**5)
    assert data["highest_jet_contraction"] < sp.Rational(3, 25)
    assert data["quadratic_Frechet_cap"] > 0
    assert data["higher_Frechet_cap"] > 0


@pytest.mark.parametrize("order,derivatives", [(2, 1), (2, 5), (3, 4), (5, 3)])
def test_independent_multinomial_allocation(order, derivatives):
    terms = [v for v in itertools.product(range(derivatives+1), repeat=order)
             if sum(v) == derivatives]
    actual_weight = sum(factorial(derivatives)//sp.prod(factorial(j) for j in v) for v in terms)
    actual_top = sum(1 for v in terms if max(v) == derivatives)
    checked = regularity.insertion_counts(order, derivatives)
    assert checked["total_multinomial_weight"] == actual_weight
    assert checked["top_jet_terms"] == actual_top
    assert checked["lower_jet_weight"] == actual_weight-actual_top


@pytest.mark.parametrize("bad", [True, 0.1, sp.Float("0.1"), "nan", "oo"])
def test_inexact_caps_rejected(bad):
    with pytest.raises((TypeError, ValueError)):
        regularity.highest_jet_bound(bad, 3, 1, 1)


@pytest.mark.parametrize("data", [(-1, 3, 1, 1), (1, 0, 1, 1), (1, 1, 1, 1),
                                  (0, 1, 1, 1), (0, 1, -1, 1)])
def test_unproved_smoothness_domain_rejected(data):
    with pytest.raises(ValueError):
        regularity.highest_jet_bound(*data)


@pytest.mark.parametrize("data", [(1, 1), (2, 0), (sp.Rational(5, 2), 2)])
def test_invalid_derivative_allocations_rejected(data):
    with pytest.raises(ValueError):
        regularity.insertion_counts(*data)


def test_regularity_auxiliary_identities():
    assert set(regularity.identities().values()) == {0}


def test_operator_norm_continuity_has_finite_actual_cap():
    assert regularity.frechet_continuity_bound(sp.Rational(1, 10**5), 3) < 10
    assert regularity.frechet_continuity_bound(0, 3) == 9
    with pytest.raises(ValueError):
        regularity.frechet_continuity_bound(1, 3)
