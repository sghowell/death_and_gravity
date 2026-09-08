"""Independent tensor contractions, Schur signs, null equations, and controls."""

from itertools import product

import pytest
import sympy as sp
from p8_affine import connection as old
from p8_affine_kinetic import vector as previous
from p8_affine_traces import geometry as g


def _zero(value):
    if isinstance(value, sp.MatrixBase):
        return all(sp.factor(entry) == 0 for entry in value)
    return sp.factor(value) == 0


def test_all_complete_scientific_identities():
    result = g.checks()
    assert len(result) == 29
    assert sum(len(value) if isinstance(value, sp.MatrixBase) else 1 for value in result.values()) == 1255
    assert all(_zero(value) for value in result.values())


def test_dense_unrestricted_two_trace_contraction_and_projective_shift():
    entries = {(a, b, c): sp.Rational((11*a+7*b+3*c) % 19-9, 13)
               for a, b, c in product(range(4), repeat=3)}
    expected_v = sp.Matrix([sum(entries[a, mu, a]-entries[a, a, mu]/4 for a in range(4))
                            for mu in range(4)])
    expected_u = sp.Matrix([sum(old.SIGN[mu]*old.SIGN[a]*entries[mu, a, a]
                                -entries[a, a, mu]/4 for a in range(4)) for mu in range(4)])
    k = sp.Matrix([entries[index] for index in old.INDICES])
    assert g.traces()["V"]*k == expected_v
    assert g.traces()["U"]*k == expected_u
    gauge = old.quadratic()["gauge"]*sp.Matrix([2, -3, 5, -7])
    assert g.traces()["full"]*(k+gauge) == expected_v.col_join(expected_u)
    assert any(entries[a, b, c] != entries[a, c, b] for a, b, c in old.INDICES)


def test_two_traces_are_covectors_under_an_exact_nontrivial_Lorentz_boost():
    boost = sp.eye(4)
    boost[:2, :2] = sp.Matrix([[5, 4], [4, 5]])/3
    metric = sp.diag(-1, 1, 1, 1)
    assert boost.T*metric*boost == metric
    inverse = boost.inv()
    distortion_transform = sp.kronecker_product(boost, inverse.T, inverse.T)
    covector_transform = sp.diag(inverse.T, inverse.T)
    assert g.traces()["full"]*distortion_transform == covector_transform*g.traces()["full"]


def test_center_eight_Schur_from_independent_full_dense_inverse():
    q = old.quotient()
    m = q["hessian"].subs(old.P, sp.Rational(1, 2))
    n = g.traces()["full"]*q["embedding"]
    d = n*m.inv(method="DM")*n.T
    expected = sp.Rational(3, 8)*sp.kronecker_product(sp.Matrix([[1, -7], [-7, 1]]),
                                                   sp.diag(1, -1, -1, -1))
    assert d == expected == g.schur()["D_clock"]
    assert d.rank() == 8
    assert d.det() == sp.Rational(531441, 256)
    assert g.traces()["source_metric"] == -g.ETA


def test_general_p_Schur_formula_is_not_an_assumed_clock_extension():
    p = old.P
    d = g.schur()["D"]
    assert sp.factor(d[4, 4]-3*(p-1)**2/2) == 0
    assert sp.factor(d[0, 4]-3*(2*p**3-2*p-1)/(4*p)) == 0
    assert sp.factor(d[1, 5]+(16*p**2-40*p-5)/(32*p**2)) == 0
    assert sp.factor(d[0, 0]-d[4, 4]) != 0


@pytest.mark.parametrize(("a", "b"), [(1, 0), (0, 1), (1, 1), (2, -3),
                                        (7+4*sp.sqrt(3), 1), (7-4*sp.sqrt(3), 1)])
def test_rank_four_is_not_inferred_from_nonzero_Schur(a, b):
    data = g.rank_one()
    n = data["N"].subs({g.A: a, g.B: b})
    assert n.rank() == 4
    params = g.require_parameters(a, b)
    d = data["D"].subs({g.A: a, g.B: b}).applyfunc(sp.simplify)
    assert d == params["gamma"]*g.ETA
    if params["gamma"] == 0:
        assert d.rank() == 0


def test_full_U_Hessian_projection_and_total_CD_coefficient_jets():
    data = g.covariant()
    assert _zero(data["reconstruction_residual"])
    assert data["clock_values"] == {"a": 0, "b": 5/(2*g.H), "c": 0, "d": 1/g.H}
    assert data["clock_X_derivatives"]["a"] == -9*g.HP/(8*g.H**2)
    assert data["clock_X_derivatives"]["c"] == 0
    assert all(value == 0 for value in data["coefficient_residuals"].values())


def test_U_rolling_lapse_coefficient_keeps_all_background_terms():
    data = g.rolling()
    s = data["symbols"]
    expected = 3*s["n_dot"]/(2*g.H)-9*g.HP*s["n"]/(4*g.H**2)
    assert sp.factor(data["U_first_variation_before_clock"][0]-expected) == 0
    total = data["U_first_variation"]+previous.adm_linearization()["first_variation_CD"]
    assert sp.factor(total[0]+15*s["Hubble"]*s["n"]/(4*g.H)) == 0
    assert total[1:, :] == sp.zeros(3, 1)
    # Ustar=-Vstar is false away from the bounce even though the derivative pieces cancel.
    assert sp.factor(total[0].subs({s["Hubble"]: 1, s["n"]: 1})) != 0
    forbidden = {s["zeta"], s["zeta_dot"], *s["shift"], *s["shift_dot"],
                 *s["shift_gradient"], *s["zeta_gradient"]}
    assert not any(entry.free_symbols & forbidden for entry in data["U_first_variation"])


def test_whole_product_gradient_shift_for_arbitrary_trace_combination():
    u, position = sp.symbols("u position", real=True)
    h = (1+u**2)**3
    hubble = 4*u/(1+u**2)
    alpha = 3*(g.A-g.B)/(2*h)
    temporal = -3*(g.A+9*g.B)*hubble/(8*h)
    assert sp.factor(-sp.diff(alpha, u)-temporal-3*(7*g.A+3*g.B)*hubble/(8*h)) == 0
    lapse = sp.Function("n")(u, position)
    scalar = alpha*lapse
    assert sp.diff(scalar, u, position)-sp.diff(scalar, position, u) == 0
    wrong_curl = sp.diff(alpha*sp.diff(lapse, position), u)-sp.diff(alpha*sp.diff(lapse, u), position)
    assert sp.factor(wrong_curl-sp.diff(alpha, u)*sp.diff(lapse, position)) == 0


def test_general_multiplier_action_and_each_full_connection_Euler_equation():
    data = g.action_identity()
    assert data["full_Euler_residual"].shape == (64, 1)
    assert data["quotient_Euler_residual"].shape == (60, 1)
    assert all(_zero(value) for name, value in data.items() if name.endswith("residual"))
    # Ordinary gamma!=0 elimination must have the PLUS Schur sign.
    r = g.rank_one()
    delta, multiplier = data["Delta"], data["multiplier"]
    inverse = g.ETA/r["gamma"]
    reduced = data["reduced_multiplier_action"].subs(dict(zip(multiplier, inverse*delta, strict=True)))
    assert sp.factor(reduced-(delta.T*inverse*delta)[0]/2) == 0


@pytest.mark.parametrize("root", [7-4*sp.sqrt(3), 7+4*sp.sqrt(3)])
def test_zero_Schur_radical_and_all_Euler_substitution(root):
    r = g.rank_one()
    substitutions = {g.A: root, g.B: 1}
    m = r["M"]
    n = r["N"].subs(substitutions)
    w = r["W"].subs(substitutions)
    assert w.rank() == 4
    assert (n*w).applyfunc(sp.simplify) == sp.zeros(4)
    assert (w.T*m*w).applyfunc(sp.simplify) == sp.zeros(4)
    arbitrary_euler = sp.Matrix([2, -3, 5, -7])
    y = -w*arbitrary_euler
    assert (m*y+n.T*arbitrary_euler).applyfunc(sp.simplify) == sp.zeros(60, 1)
    assert (n*y).applyfunc(sp.simplify) == sp.zeros(4, 1)
    assert sp.simplify((y.T*m*y)[0]) == 0
    assert sp.simplify(root-1).is_nonzero is True


def test_null_Schur_is_only_quadratic_rolling_not_nonlinear_open_tube():
    a = 7+4*sp.sqrt(3)
    combine = sp.eye(4).row_join(sp.eye(4))
    combine[:, :4] *= a
    d_off = combine*g.schur()["D"].subs(old.P, sp.Rational(51, 100))*combine.T
    assert sp.simplify(d_off[0, 0]).is_nonzero is True
    null = g.null_schur()
    assert all(_zero(value) for name, value in null.items() if name.endswith("remainder"))
    assert (null["kernel_dimension"], null["kernel_radical_dimension"],
            null["remaining_nondegenerate_kernel_quotient_dimension"]) == (56, 4, 52)
    assert null["nonlinear_open_tube_inverse_claim"] is False


def test_null_action_is_the_actual_Maxwell_curl_of_the_rolling_lapse_source():
    data = g.rolling()
    null = g.null_schur()
    s = data["symbols"]
    curl = (-(g.A-g.B)*s["n_dot"]/g.H
            +3*(11*g.A-g.B)*s["Hubble"]*s["n"]/(8*g.H))
    expected = null["kinetic_coefficient"]*null["q_spatial"]*curl**2/2
    assert sp.factor(null["Maxwell_scalar_density_divided_by_a_cubed"]-expected) == 0
    assert sp.factor(null["lapse_dot_squared_coefficient"]
                     -null["kinetic_coefficient"]*null["q_spatial"]*(g.A-g.B)**2/(2*g.H**2)) == 0


@pytest.mark.parametrize("kinetic", [5, -5])
def test_nonunit_normalization_accepts_both_nonzero_curl_signs(kinetic):
    data = g.units(1, 0, mass_squared=3, time_scale=2, kinetic_coefficient=kinetic)
    assert data["lambda_normalized"] == sp.Rational(kinetic, 12)
    assert data["D_physical_factor"] == sp.Rational(1, 3)
    assert data["trace_normalized_factor"] == 2


@pytest.mark.parametrize(("a", "b", "kwargs"), [
    (0, 0, {}), (True, 1, {}), (0.5, 1, {}), ("1", 0, {}),
    (sp.oo, 1, {}), (sp.nan, 1, {}), (sp.I, 0, {}), (sp.Symbol("x"), 0, {}),
    (1, 0, {"kinetic_coefficient": 0}), (1, 0, {"kinetic_coefficient": 1.0}),
    (1, 0, {"mass_squared": -1}), (1, 0, {"time_scale": 0}),
    (1, 0, {"time_scale": "2"}),
])
def test_invalid_exact_parameters_are_rejected(a, b, kwargs):
    with pytest.raises((TypeError, ValueError)):
        g.require_parameters(a, b, **kwargs)
