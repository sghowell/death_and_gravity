"""Two-frequency constants, exact proper clock and sampler-domain checks."""

import pytest
import sympy as sp
from p8a_see_qsei import sampling


def test_exact_spectral_and_clock_identities():
    assert all(value == 0 for value in sampling.spectral_identities().values())
    assert all(value == 0 for value in sampling.clock_identities().values())


def test_all_error_pieces_are_positive_and_accounted_for():
    data = sampling.calibration()
    pieces = data["spectral_error_roots"]
    assert set(pieces) == {"forward", "local", "real_history_Parseval", "history_remainder", "infrared"}
    assert all(value > 0 for value in pieces.values())
    assert data["total_root"] == data["auxiliary_root"]+sum(pieces.values())
    assert data["total_root"] > data["auxiliary_root"]


def test_exact_root_coefficient_and_coarsening():
    data = sampling.calibration()
    root = sp.Rational(6400000017560000000992361942917377673,
                       6400000000000000000000000000000000000)
    assert data["total_root"] == root
    assert 1 < root < 1+sp.Rational(1, 10**8)
    assert data["derived_coefficient"] == root**2
    assert data["derived_coefficient"] < data["rounded_coefficient"] == 2
    assert all(value > 0 for value in data["strict_margins"].values())


def test_named_bound_magnitude_has_exact_physical_scaling():
    assert sampling.absolute_bound(0) == 0
    assert sampling.absolute_bound(3, hbar=4) == sp.Rational(3, 2)/sp.pi**2
    assert sampling.absolute_bound(7, hbar=2) == 2*sampling.absolute_bound(7)
    assert sampling.absolute_bound(14) == 2*sampling.absolute_bound(7)


@pytest.mark.parametrize("norm,hbar", [(True, 1), (0.1, 1), (-1, 1), (1, sp.oo), (1, -1), (1, False)])
def test_bound_rejects_inexact_or_wrong_sign_inputs(norm, hbar):
    with pytest.raises((TypeError, ValueError)):
        sampling.absolute_bound(norm, hbar=hbar)


def test_constant_scale_clock_exposes_missing_measure_factor():
    x = sp.Symbol("x", real=True)
    scale = sp.Integer(3)
    psi = sp.Function("psi")
    f = scale**-sp.Rational(3, 2)*psi(scale*x)
    # F''/sqrt(a)=psi_ss, not F'' itself: L2 changes by dx=ds/a.
    target = sp.diff(psi(x), x, 2).subs(x, scale*x)
    assert sp.simplify(sp.diff(f, x, 2)/sp.sqrt(scale)-target) == 0
    assert sp.simplify(sp.diff(f, x, 2)-target) != 0


def test_auxiliary_positive_form_needs_cross_term():
    first, weighted = sp.symbols("first weighted", real=True)
    correct = first**2+2*weighted**2+sp.Rational(8, 3)*first*weighted
    omitted = first**2+2*weighted**2
    assert correct.subs({first: 1, weighted: 1}) > omitted.subs({first: 1, weighted: 1})
    assert sp.Rational(3, 2)**2 >= 2
