from itertools import combinations

import pytest
import sympy as sp
from p8_m1_physical import lagrangian, regressions, vertices
from p8_m1_physical.vertices import Leg, tensor_basis


def test_fixture_all_external_and_internal_gamma_modes_are_admitted():
    for waves in (regressions.TRI, regressions.QUAD):
        assert tuple(map(sum, zip(*waves))) == (0, 0, 0)
        for n in range(1, len(waves)):
            for indices in combinations(range(len(waves)), n):
                transfer = tuple(sum(waves[i][j] for i in indices) for j in range(3))
                assert sum(k*k for k in transfer) > 6


@pytest.mark.parametrize("chart,point", [("gamma", 0), ("unitary", sp.Rational(1, 3))])
@pytest.mark.parametrize("kinds", [("s", "m", "P"), ("s", "m", "t", "pi"),
                                    ("p", "P", "p", "P"), ("m", "m", "m", "m")])
def test_constraints_reality_parity_permutation(chart, point, kinds):
    legs = regressions.fixture(kinds)
    result = vertices.hamiltonian_kernel(legs, point, chart)
    assert all(result["constraint_checks"].values())
    assert sp.im(result["kernel"]) == 0
    reversed_legs = tuple(reversed(legs))
    parity = tuple(Leg(tuple(-k for k in leg.wave), leg.kind, leg.polarization) for leg in legs)
    assert vertices.hamiltonian_kernel(reversed_legs, point, chart)["kernel"] == result["kernel"]
    assert vertices.hamiltonian_kernel(parity, point, chart)["kernel"] == result["kernel"]


def test_cubic_and_quartic_literal_bounce_regressions():
    for kinds, expected in {
        ("s", "m", "m"): sp.Rational(46767803, 239800),
        ("s", "s", "m", "m"): sp.Rational(399455066903941, 86184179950),
        ("m", "m", "m", "m"): -sp.Rational(2775087, 43600),
    }.items():
        assert vertices.hamiltonian_kernel(regressions.fixture(kinds), 0)["kernel"] == expected


def test_quartic_requires_full_scalar_matrix_contact():
    legs = regressions.fixture(("m_dot",)*4)
    result = lagrangian.kernel(legs, 0, "gamma")
    assert result["kernel"] == sp.Rational(15739581640021, 33058771200000)
    assert result["minus_H"] == -sp.Rational(3269548239, 12940480000)
    assert all(item["scalar_off_diagonal"] != 0 for item in result["Legendre_corrections"])
    assert result["kernel"] == result["minus_H"]+sum(item["scalar"]+item["tensor"]
                                                     for item in result["Legendre_corrections"])
    assert lagrangian.kernel(tuple(reversed(legs)), 0, "gamma")["kernel"] == result["kernel"]


def test_both_tensor_polarizations_and_internal_basis_invariance(monkeypatch):
    legs = regressions.fixture(("s", "m_dot", "s", "m_dot"))
    result = lagrangian.kernel(legs, 0, "gamma")
    assert any(all(value != 0 for value in item["tensor_polarizations"])
               for item in result["Legendre_corrections"])

    def rotated_basis(wave):
        E, F = map(sp.Matrix, tensor_basis(wave))
        ratio = sp.trace(F*F)/sp.trace(E*E)
        return tuple(tuple(tuple(row) for row in matrix.tolist()) for matrix in (E+F, ratio*E-F))

    try:
        with monkeypatch.context() as patch:
            patch.setattr(lagrangian, "tensor_basis", rotated_basis)
            lagrangian.kernel.cache_clear()
            rotated = lagrangian.kernel(legs, 0, "gamma")
            assert rotated["kernel"] == result["kernel"]
            assert [item["tensor"] for item in rotated["Legendre_corrections"]] == [
                item["tensor"] for item in result["Legendre_corrections"]]
    finally:
        lagrangian.kernel.cache_clear()


def test_exceptional_and_invalid_inputs_are_rejected():
    with pytest.raises(ValueError):
        vertices.hamiltonian_kernel((Leg((1, 0, 0), "s"), Leg((1, 0, 0), "m")), 0)
    with pytest.raises(ValueError):
        vertices.hamiltonian_kernel(tuple(Leg(k, "m") for k in ((1, 0, 0), (-1, 0, 0),
                                                               (0, 1, 0), (0, -1, 0))), 0)
    with pytest.raises(ValueError):
        vertices.hamiltonian_kernel(regressions.fixture(("s", "m", "P")), 2)
    with pytest.raises(ValueError):
        vertices.hamiltonian_kernel((Leg((0, 0, 0), "s"),), 0, "gamma")
    with pytest.raises(ValueError):
        lagrangian.kernel(regressions.fixture(("p", "m", "P")), 0)
    bad = tuple(tuple(row) for row in sp.eye(3).tolist())
    with pytest.raises(ValueError):
        vertices.hamiltonian_kernel((Leg((2, 2, 0), "t", bad), Leg((-2, -2, 0), "t", bad)), 0)
