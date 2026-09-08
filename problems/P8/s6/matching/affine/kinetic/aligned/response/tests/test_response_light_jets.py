"""Actual ODE source, continuous coefficient majorants and convolution controls."""
import pytest
import sympy as sp
from p8_aligned_response import leading as l


def test_actual_ODE_second_variation_and_preparation():
    assert all(value == 0 for value in l.checks().values())
    assert l.source()["Q_second"] != 0
    assert l.source()["expected"].subs({l.N: 1, l.K: 1, l.u: 0}) == -2


def test_exact_jet_constants_and_source_envelope():
    data = l.jet_constants()
    assert data["mixed_source_jets"] == [2, 10, 65, 514]
    assert data["squared_source_jets"] == [sp.Rational(483, 16), sp.Rational(6429, 32),
                                          sp.Rational(12921, 8), sp.Rational(132027, 8)]
    assert l.source_envelope(1, 1)["source_jet_E"] == sp.Rational(136139, 8)
    assert l.source_envelope(0, 1)["source_jet_E"] == 0
    assert l.source_envelope(1, 1)["bound_on_full_nonlinear_source_or_on_shell_evolution"] is False


@pytest.mark.parametrize("point", (-sp.Rational(1, 2), -sp.Rational(1, 4), 0, sp.Rational(1, 4), sp.Rational(1, 2)))
def test_independent_coefficient_derivative_bounds(point):
    data, limits = l.source(), l.jet_constants()
    for coefficient, key in ((data["mixed_coefficient"], "A_coefficient_jets"),
                             (data["squared_lapse_coefficient"], "C_coefficient_jets")):
        for j in range(4):
            assert abs(sp.diff(coefficient, l.u, j).subs(l.u, point)) <= limits[key][j]


def test_fourier_source_uses_output_convolution_not_input_momentum():
    incoming = {-1: sp.Rational(1, 2), 1: sp.Rational(1, 2)}
    product = {}
    for k, a in incoming.items():
        for j, b in incoming.items():
            product[k+j] = product.get(k+j, 0)+a*b
    assert product == {-2: sp.Rational(1, 4), 0: sp.Rational(1, 2), 2: sp.Rational(1, 4)}
    assert 1 not in product
    assert sum(abs(value) for value in product.values()) <= sum(abs(value) for value in incoming.values())**2


@pytest.mark.parametrize("values", ((-1, 1), (1, -1), (0.1, 1), (True, 1), (1, sp.oo)))
def test_invalid_light_jet_envelopes_are_rejected(values):
    with pytest.raises((TypeError, ValueError)):
        l.source_envelope(*values)
