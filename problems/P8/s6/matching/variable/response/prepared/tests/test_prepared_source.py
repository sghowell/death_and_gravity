from fractions import Fraction

import pytest
import sympy as sp
from p8_variable_prepared import independent, source


def test_exact_pole_cancelled_full_loading_identity():
    checks = source.checks()
    assert checks["full_pole_cancelled_commutator_loading"] == 0
    assert checks["fixed_cutoff_transition_left"] == checks["fixed_cutoff_transition_right"] == 0
    assert checks["source_prefactor_margin"] > 0


def test_independent_source_jet_calibration_and_norm():
    separate = independent.source_constants()
    assert source.coefficient_bounds() == {key: sp.Rational(value) for key, value in separate["coefficients"].items()}
    switch = source.switch_bounds()
    for name in ("bump_jets", "reciprocal_jets"):
        assert switch[name] == [sp.Rational(value) for value in separate[name]]
    assert switch["zeta_jets"] == {key: sp.Rational(value) for key, value in separate["zeta_jets"].items()}
    assert source.source_sup_bound() == separate["source_sup"] == 321137485366608000
    assert switch["bump_jets"] == [1, 2, 36, 1584, 129600]


def test_delta_dependent_source_retuning_has_fixed_support_and_norm():
    delta = sp.Rational(1, 10**12)
    assert source.source_difference_bound(delta, 2) == 800*delta*source.source_sup_bound()
    assert source.source_sup_bound(2) == 2*source.source_sup_bound()
    assert source.switch_bounds()["support_length"] == sp.Rational(1, 100)


def test_fourth_cutoff_derivative_cannot_be_omitted():
    u = sp.Symbol("u")
    values = [sp.Function(name)(u) for name in ("kg", "kf", "gg", "P", "g", "f")]
    data = source.loading_coefficients(u, *values, sp.Symbol("K"))
    kg, kf, _, inv, _, f = values
    assert sp.expand(data["coefficients"][4]-kg*inv*kf*f) == 0
    assert data["coefficients"][4] != 0


def test_generic_cutoff_interface_does_not_assume_calibrated_shape():
    jets = {1: 1, 2: 2, 3: 3, 4: 4}
    expected = sum(source.coefficient_bounds()[key]*value for key, value in jets.items())
    assert source.source_sup_bound(1, jets) == expected


@pytest.mark.parametrize("call", [
    lambda: source.source_sup_bound(True),
    lambda: source.source_sup_bound(1, {1: 1, 2: 2, 3: 3}),
    lambda: source.source_sup_bound(1, {1: 1, 2: -2, 3: 3, 4: 4}),
    lambda: source.source_sup_bound(1, {1: 1, 2: 2, 3: 3, 4: 4.0}),
    lambda: source.source_difference_bound(Fraction(1, 100), 1),
])
def test_invalid_source_norms_and_missing_jets_are_rejected(call):
    with pytest.raises((TypeError, ValueError)):
        call()
