"""Exact macroscopic parameter bounds and explicit source/prescription gates."""

import pytest
import sympy as sp
from p8a_maxwell_focusing import envelope, focusing


def test_all_scaling_and_physical_dictionary_identities():
    assert all(sp.simplify(value) == 0 for value in envelope.identities().values())
    assert all(sp.simplify(value) == 0 for value in focusing.identities().values())
    assert all(value > 0 for value in envelope.radical_margins().values())


def test_demonstration_costs_and_margins_are_exact():
    data = focusing.calibration()
    base = data["zero_source"]
    assert base["envelope"]["past_cost"] == sp.Rational(54968328329, 4375)
    assert base["envelope"]["future_cost"] == sp.Rational(363631, 240)
    assert base["envelope"]["total_cost"] < 13000000
    assert base["strict_focusing_margin"] == sp.Rational(2398243151869, 3000000000000)
    assert data["source_sigma_at_most_one"]["strict_focusing_margin"] > sp.Rational(1, 3)


def test_large_source_penalty_prevents_false_focusing():
    data = focusing.theorem_constants(sp.Rational(1, 100), sp.Rational(3, 2),
                                      (2, 4, 16, 96), (2, 4, 16, 96),
                                      beta_m=1, delta=sp.Rational(1, 10**8), sigma=100)
    assert data["sufficient_focusing_test"] is False
    assert data["index_plus_K_times_tau_upper"] > 0


def test_signed_beta_is_visible_but_the_envelope_uses_its_magnitude():
    positive = envelope.coefficients((2, 4, 16, 96), (2, 4, 16, 96), sp.Rational(1, 100), beta_m=1)
    negative = envelope.coefficients((2, 4, 16, 96), (2, 4, 16, 96), sp.Rational(1, 100), beta_m=-1)
    assert positive["beta_M"] == -negative["beta_M"] == 1
    assert positive["total_cost"] == negative["total_cost"]
    with pytest.raises(TypeError):
        envelope.coefficients((2, 4, 16, 96), (2, 4, 16, 96), 1)


@pytest.mark.parametrize("bad", [True, 0.1, sp.oo, sp.nan, sp.Symbol("x")])
def test_delta_input_is_exact_finite_and_resolved(bad):
    with pytest.raises((TypeError, ValueError)):
        focusing.theorem_constants(1, 1, (2, 4, 16, 96), (2, 4, 16, 96),
                                    beta_m=0, delta=bad, sigma=0)


def test_incompatible_caps_are_rejected_not_used_as_a_vacuous_witness():
    with pytest.raises(ValueError, match="contradicts"):
        focusing.theorem_constants(1, 2, (1, 4, 16, 96), (2, 4, 16, 96),
                                    beta_m=0, delta=0, sigma=0)
    with pytest.raises(TypeError, match="four exact"):
        envelope.coefficients((1, 2, 3), (1, 2, 3, 4), 1, beta_m=0)


def test_physical_enclosure_keeps_all_source_signs():
    positive = focusing.physical_enclosure(kappa=2, hbar=3, tau=5,
                                            cosmological_constant=7, other_eed_lower=-11)
    negative = focusing.physical_enclosure(kappa=2, hbar=3, tau=5,
                                            cosmological_constant=-7, other_eed_lower=11)
    assert positive["sigma"] == 25*29
    assert negative["sigma"] == 0
    assert positive["rational_delta_upper"] == sp.Rational(1, 300)
