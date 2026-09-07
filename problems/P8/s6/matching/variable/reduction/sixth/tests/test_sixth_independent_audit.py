"""Separately authored coordinate, ADM, root-series and cancellation checks."""
import pytest
import sympy as sp
from p8_sixth_reduction import bridges, literal, tensor


def test_all_forty_independent_and_continuous_identities():
    values = bridges.checks()
    assert len(values) == 40
    assert set(values.values()) == {0}


def test_full_coordinate_geometry_includes_lapse_and_trace():
    data = literal.coordinate_schouten()
    # At q'=q''=ell=0, the time component still contains -2 ell'.
    h0 = data["h"][0].subs({literal.Q1: 0, literal.Q2: 0, literal.L: 0})
    assert h0 == -2*literal.A*literal.L1
    assert data["mixed_correction_over_A"].shape == (4, 4)
    assert data["scalar"].has(sp.diff(literal.R, literal.TIME, 2))


@pytest.mark.parametrize("n", (1, 2, 3))
def test_independent_normal_coefficient(n):
    independent = literal.normal_and_euler()[f"B{n}"].subs(bridges.literal_to_primary())
    assert sp.factor(independent-tensor.derive()[f"B{n}"]) == 0


@pytest.mark.parametrize("n", (2, 4, 6))
def test_physical_scale_dictionary_from_independent_center(n):
    point, primary = literal.independent_center(), tensor.center()
    expected = tensor.M**2*tensor.TAU**(n-2)*point["E"][n].subs(point["c"], tensor.C)
    assert sp.factor(primary["sixth_E_physical"][n]-expected) == 0
    assert sp.diff(expected, tensor.TAU)*tensor.TAU == (n-2)*expected


def test_compact_probe_jet_and_strict_continuous_individual_floor():
    point = literal.independent_center()
    delta = sp.Symbol("delta", positive=True)
    e2 = sp.Poly(sp.expand(point["E"][2].subs(point["c"], 2+delta)), delta)
    assert e2.nth(0) == -2
    assert e2.nth(1) == sp.Rational(7, 4)
    assert all(e2.nth(n) < 0 for n in (2, 3, 4))
    assert 2-e2.nth(1)/100 == sp.Rational(793, 400)
    # Any smooth compact q agreeing with u^2/2 near zero has this jet;
    # its q^(4) and q^(6) entries vanish at the center, not globally.
    u = sp.Symbol("u", real=True)
    q = u*u/2
    values = {n: sp.diff(q, u, n).subs(u, 0) for n in range(7)}
    assert sum(point["E"][n]*values[n] for n in range(7)) == point["E"][2]


def test_individual_nonzero_term_does_not_survive_the_required_cancellation():
    point = literal.independent_center()
    c = point["c"]
    e4_low = (9*c*c-22*c+24)/8
    assert sp.limit(point["E"][2], c, 2, dir="+") == -2
    assert sp.limit(e4_low, c, 2, dir="+") == 2
    combined = sp.factor(point["E"][2]+e4_low)
    assert sp.limit(combined, c, 2, dir="+") == 0
    assert sp.limit(combined/(c-2), c, 2, dir="+") == sp.Rational(7, 2)
    assert combined != 0  # cancellation of the limit is not an all-c identity
