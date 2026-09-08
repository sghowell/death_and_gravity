"""Independent spatial Legendre maps and omitted-Gauss-term controls."""
import pytest
import sympy as sp
from p8_affine_nonlinear import maxwell, remaining


def test_complete_vector_and_remaining_velocity_identities():
    for module in (maxwell, remaining):
        for value in module.checks().values():
            assert all(item == 0 for item in value) if isinstance(value, sp.MatrixBase) else value == 0


@pytest.mark.parametrize("G", (sp.eye(3), sp.Matrix([[2, 1, 0], [1, 3, 1], [0, 1, 2]])))
@pytest.mark.parametrize("N,kappa", ((sp.S.One, sp.S.One), (sp.Rational(3, 2), sp.Rational(1, 7))))
def test_independent_dense_spatial_Maxwell_Legendre_solve(G, N, kappa):
    assert all(G[:k, :k].det() > 0 for k in range(1, 4))
    data = maxwell.legendre()
    E, P, DW0 = [data[name] for name in ("E", "P", "DW0")]
    L = kappa*(E.T*G*E)[0]/(2*N)
    solution = sp.solve([sp.diff(L, E[i])-P[i] for i in range(3)], list(E))
    actual = sp.expand(((P.T*(E+DW0))[0]-L).subs(solution))
    expected = N*(P.T*G.inv()*P)[0]/(2*kappa)+(P.T*DW0)[0]
    assert sp.expand(actual-expected) == 0
    assert not (actual-(P.T*DW0)[0]).has(*DW0)
    assert sp.hessian(L, list(E)).det() == (kappa/N)**3*G.det()


def test_negative_curl_sign_is_not_in_the_positive_velocity_branch():
    data = maxwell.legendre()
    E = data["E"]
    L = -(E.T*E)[0]/2
    assert sp.hessian(L, list(E)) == -sp.eye(3)


def test_spatial_generator_includes_temporal_Gauss_contribution():
    data = maxwell.spatial_generator()
    xi, W = [sp.Matrix(sp.symbols(prefix+"0:3", real=True)) for prefix in ("xi", "W")]
    divP = sp.Symbol("div_Pi", real=True)
    omitted = sp.expand(data["target"]+(xi.T*W)[0]*divP)
    assert sp.expand(omitted-data["target"]-(xi.T*W)[0]*divP) == 0
    assert sp.expand(omitted-data["target"]) != 0


def test_all_five_shear_and_matter_velocity_pivots_are_retained():
    shear, matter = remaining.shear(), remaining.matter()
    assert sp.hessian(shear["L"], list(shear["v"])) == 2*shear["U"]*shear["minus_B"]*sp.eye(5)
    assert sp.diff(matter["H"], matter["p"], 2) == matter["N"]/matter["U"]
