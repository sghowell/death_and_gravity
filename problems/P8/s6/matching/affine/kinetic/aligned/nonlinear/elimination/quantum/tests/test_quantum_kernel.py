"""Independent full determinant and physical longitudinal controls."""
import pytest
import sympy as sp
from p8_aligned_quantum import kernel


def test_full_determinant_and_temporal_constraint_identities():
    assert all(value == 0 for value in kernel.checks().values())


@pytest.mark.parametrize("a,b", ((1, 1), (sp.Rational(9, 10), sp.Rational(21, 20)), (2, 3)))
def test_independent_dense_momentum_rank_one_determinant(a, b):
    # No special alignment of spatial momentum is used in this control.
    p = sp.Matrix([2, 1, 3, -1])
    mass = sp.diag(a, b, b, b)
    K = (p.T*p)[0]*sp.eye(4)-p*p.T+mass
    expected = (4+11+b)**2*(a*b+4*a+11*b)
    assert K.det() == expected
    diagonal = (p.T*p)[0]*sp.eye(4)+mass
    lemma = diagonal.det()*(1-(p.T*diagonal.inv()*p)[0])
    assert sp.factor(K.det()-lemma) == 0


def test_direct_temporal_solve_matches_the_physical_dispersion():
    d = kernel.longitudinal()
    solved = sp.solve(sp.diff(d["L"], d["t"]), d["t"])
    assert len(solved) == 1
    assert sp.factor(solved[0]-d["solution"]) == 0
    reduced = d["L"].subs(d["t"], solved[0])
    omega2 = -sp.diff(reduced, d["sigma"], 2)/sp.diff(reduced, d["v"], 2)
    assert sp.factor(omega2-d["frequency"]) == 0


def test_counting_only_transverse_modes_loses_a_momentum_dependent_factor():
    d = kernel.determinant()
    difference = sp.factor(d["determinant"]/(kernel.a*kernel.m2)
                          -(kernel.omega**2+kernel.k**2+kernel.b*kernel.m2)**2)
    assert difference != 0
    assert difference.has(kernel.omega, kernel.k)


def test_homogeneous_limit_keeps_three_massive_polarizations():
    d = kernel.determinant()
    expected = kernel.a*kernel.m2*(kernel.omega**2+kernel.b*kernel.m2)**3
    assert sp.factor(d["determinant"].subs(kernel.k, 0)-expected) == 0
