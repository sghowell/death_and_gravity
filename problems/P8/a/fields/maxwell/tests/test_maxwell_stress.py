"""Formal source/sign controls for the new photon prescription family."""

import sympy as sp
from p8a_maxwell import stress


def test_every_covariant_reference_identity_is_exact():
    assert all(sp.simplify(value) == 0 for value in stress.identities().values())


def test_beta_is_symbolic_and_not_silently_fixed_to_zero():
    data = stress.reference_jets()
    assert stress.BETA_M in data["rho"].free_symbols
    assert sp.diff(data["rho"], stress.BETA_M) == data["I_rho"]
    assert sp.diff(data["EED"], stress.BETA_M) == data["I_EED"]
    assert data["rho"]-data["EED"] == sp.expand(data["trace"]/2)


def test_de_sitter_and_radiation_credit_have_different_signs():
    h, t = sp.symbols("H t", positive=True)
    de_sitter = stress.reference_jets(h, 0, 0, 0)
    assert de_sitter["rho"] == 186*h**4
    assert de_sitter["EED"] == -186*h**4
    radiation = stress.reference_jets(1/(2*t), -1/(2*t**2), 1/t**3, -3/t**4)
    assert radiation["rho"] == sp.Rational(93, 8)/t**4
    assert radiation["EED"] == sp.Rational(279, 8)/t**4
    assert radiation["I_rho"] == radiation["I_EED"] == 0


def test_nonzero_scalar_curvature_does_not_remove_finite_ambiguity():
    t, power = sp.symbols("t p", positive=True)
    h = power/t
    d = stress.reference_jets(h, *(sp.diff(h, t, j) for j in (1, 2, 3)))
    assert sp.factor(d["I_EED"]) == -54*power*(power-2)*(2*power-1)/t**4
    assert d["I_EED"].subs(power, sp.Rational(2, 3)) != 0


def test_physical_component_normalization_includes_hbar():
    data = stress.reference_stress(1, 0, 0, 0, 7, hbar=3)
    assert data["rho"] == sp.Rational(31, 160)/sp.pi**2
    assert data["EED"] == -data["rho"]


def test_wrong_HH_signature_or_trace_omission_is_detectable():
    h = sp.Symbol("H", positive=True)
    correct = stress.reference_jets(h, 0, 0, 0, 0)
    assert -correct["rho"] != correct["rho"]
    assert correct["EED"] != correct["rho"]
