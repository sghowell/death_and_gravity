"""Actual coefficient jets and an independent coordinate commutator check."""
import sympy as sp
from p8_affine import connection
from p8_affine_ricci import geometry, source


def test_all_full_source_identities():
    for value in source.checks().values():
        assert all(entry == 0 for entry in value) if isinstance(value, sp.MatrixBase) else value == 0


def test_spatial_Hessian_and_shift_were_retained_before_their_cancellation():
    result = source.stationary_variation()
    assert len(result["hessian_jets"]) == 10
    for i in range(1, 4):
        for j in range(i, 4):
            assert not any(value.has(result["hessian_jets"][connection.H[i, j]]) for value in result["actual"])
    assert sum(value != 0 for value in result["actual"]) == 28
    assert result["background"] == sp.zeros(64, 1)


def test_full_six_function_coordinate_covariant_commutator():
    data = source.rolling()
    coordinates = data["coordinates"]
    gamma, inverse = data["Christoffel"], data["inverse_metric"]
    z = [sp.Function(f"Z{i}")(*coordinates) for i in range(6)]
    lower = geometry.two_form(z)
    raised = sp.diag(*inverse)*lower
    kappa = {index: sp.diff(raised[index[0], index[1]], coordinates[index[2]])+sum(
        gamma[index[0], e, index[2]]*raised[e, index[1]]
        -gamma[e, index[1], index[2]]*raised[index[0], e] for e in range(4))
        for index in connection.INDICES}
    actual = source.ricci_difference(kappa, gamma, data["metric"], inverse, coordinates)
    u = source.u
    expected = 8*(1+7*u**2)*sp.Matrix(z)/(1+u**2)**2
    assert (actual-expected).applyfunc(sp.factor) == sp.zeros(6, 1)


def test_actual_curved_source_has_no_surviving_lapse_time_derivative():
    data = source.rolling()
    for i in range(1, 4):
        expected = -10*source.u*sp.diff(data["lapse"], data["coordinates"][i])/(1+source.u**2)**4
        assert sp.factor(data["Cstar"][i-1]-expected) == 0
    assert data["Cstar"][3:, :] == sp.zeros(3, 1)
    assert not any(entry.has(sp.diff(data["lapse"], source.u)) for entry in data["Cstar"])
