"""Independent curvature contractions, Lorentz orbits and all64 controls."""
from itertools import product

import pytest
import sympy as sp
from p8_affine import connection
from p8_affine_ricci import geometry as g
from p8_affine_traces import geometry as traces


def test_all_literal_geometry_identities():
    for value in g.checks().values():
        assert all(entry == 0 for entry in value) if isinstance(value, sp.MatrixBase) else value == 0


def test_independent_dense_curvature_contraction_not_a_trace_curl_fit():
    values = {index: sp.Rational((11*index[0]+7*index[1]+3*index[2]) % 19-9, 13)
              for index in connection.INDICES}
    k = (2, -3, 5, -7)

    def curvature(a, b, c, d):
        return k[c]*values[a, b, d]-k[d]*values[a, b, c]

    difference = sp.Matrix(4, 4, lambda mu, nu: sum(
        curvature(a, nu, a, mu)-g.SIGNS[nu]*g.SIGNS[a]*curvature(nu, a, a, mu)
        for a in range(4)))
    expected = sp.Matrix([(difference[a, b]-difference[b, a])/2 for a, b in g.PAIRS])
    actual = g.flat_map().subs(dict(zip(g.K, k, strict=True)))*sp.Matrix([values[i] for i in connection.INDICES])
    assert actual == expected
    axial = sp.zeros(64, 1)
    axial[connection.index(0, 1, 2)], axial[connection.index(0, 2, 1)] = 1, -1
    assert traces.traces()["full"]*axial == sp.zeros(8, 1)
    visible = g.flat_map().subs(dict(zip(g.K, (1, 0, 0, 0), strict=True)))*axial
    assert visible[g.PAIRS.index((1, 2))] == -1


@pytest.mark.parametrize("covector", ((1, 0, 0, 0), (0, 0, 0, 1), (1, 0, 0, 1)))
def test_rank_six_for_all_three_nonzero_Lorentz_orbits(covector):
    n = g.flat_inverse()["N"].subs(dict(zip(g.K, covector, strict=True)))
    assert n.rank() == 6
    m = g.flat_inverse()["M"]
    assert n*m.inv(method="DM")*n.T == sp.zeros(6)


def test_exact_Lorentz_covariance_of_the_full_linear_curvature_map():
    boost = sp.eye(4)
    boost[:2, :2] = sp.Matrix([[5, 4], [4, 5]])/3
    inverse = boost.inv()
    metric = sp.diag(*g.SIGNS)
    assert boost.T*metric*boost == metric
    z = sp.symbols("z0:6", real=True)
    changed = inverse.T*g.two_form(z)*inverse
    pair_transform = sp.Matrix([changed[a, b] for a, b in g.PAIRS]).jacobian(z)
    lower_k = inverse.T*sp.Matrix(g.K)
    changed_map = g.flat_map().subs(dict(zip(g.K, lower_k, strict=True)), simultaneous=True)
    tensor_transform = sp.kronecker_product(boost, inverse.T, inverse.T)
    assert (changed_map*tensor_transform-pair_transform*g.flat_map()).applyfunc(sp.expand) == sp.zeros(6, 64)


def test_flat_nilpotent_update_has_an_exact_polynomial_inverse():
    data = g.flat_inverse()
    n = data["N"].subs(dict(zip(g.K, (2, 1, -3, 4), strict=True)))
    m = data["M"]
    coupling = sp.Rational(5, 12)
    # Overall Fourier/adjoint signs do not affect nilpotency; test both.
    inverse = m.inv(method="DM")
    for sign in (-1, 1):
        update = sign*coupling*n.T*g.PAIRING*n
        candidate = inverse-inverse*update*inverse
        assert (m+update)*candidate == sp.eye(60)
        assert (inverse*update)**2 == sp.zeros(60)


def test_FLRW_Riemann_tensor_has_the_required_symmetries():
    result = g.commutator()
    tensor = result["Riemann"]
    for a, b, c, d in product(range(4), repeat=4):
        assert tensor[a, b, c, d] == -tensor[a, b, d, c]
        assert g.SIGNS[a]*tensor[a, b, c, d] == -g.SIGNS[b]*tensor[b, a, c, d]
        assert g.SIGNS[a]*tensor[a, b, c, d] == g.SIGNS[c]*tensor[c, d, a, b]
    assert (result["endomorphism"]-result["scalar_curvature"]*sp.eye(6)/3).applyfunc(sp.factor) == sp.zeros(6)


@pytest.mark.parametrize("p", (sp.Rational(49, 100), sp.Rational(51, 100)))
def test_off_clock_null_identity_fails_but_does_not_imply_a_new_time_mode(p):
    data = g.off_clock_symbol()
    response = data["response"].subs({connection.P: p, data["spatial_k"]: 1})
    assert response.rank() == 4
    assert not response.has(data["omega"])
    assert data["response"].subs(data["spatial_k"], 0) == sp.zeros(6)
    assert data["additional_propagating_mode_or_ghost_claim"] is False
