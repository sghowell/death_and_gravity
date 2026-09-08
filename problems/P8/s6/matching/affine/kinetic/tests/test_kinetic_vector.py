"""Independent Schur, complete ADM pullback, domain, and omission controls."""

from itertools import product

import pytest
import sympy as sp
from p8_affine import connection as old
from p8_affine_kinetic import vector as v


def _zero(value):
    if isinstance(value, sp.MatrixBase):
        return all(entry == 0 for entry in value)
    return value == 0


def test_every_exact_residual():
    residuals = v.checks()
    assert len(residuals) == 30
    assert all(_zero(value) for value in residuals.values())


def test_vector_contraction_on_an_unrestricted_dense_connection():
    k = {(a, b, c): sp.Rational((13*a+7*b+5*c) % 17-8, 9)
         for a, b, c in product(range(4), repeat=3)}
    expected = sp.Matrix([sum(k[a, mu, a]-k[a, a, mu]/4 for a in range(4))
                          for mu in range(4)])
    assert v.vector_map()["full"]*sp.Matrix([k[key] for key in old.INDICES]) == expected
    assert v.vector_map()["quotient"].rank() == 4
    assert _zero(v.vector_map()["projective_residual"])


@pytest.mark.parametrize("p", [sp.Rational(1, 2), sp.sqrt(sp.Rational(9, 40)),
                                sp.sqrt(sp.Rational(11, 40))])
def test_schur_matrix_from_a_separate_dense_numeric_inverse(p):
    q = old.quotient()
    m = q["hessian"].subs(old.P, p)
    n = v.vector_map()["quotient"]
    independent = n*m.inv(method="DM")*n.T
    assert all(sp.simplify(value) == 0 for value in independent-v.schur()["D"].subs(v.P, p))


def test_fixed_vector_action_and_all_56_complement_directions():
    data = v.schur()
    p = sp.Rational(1, 2)
    m = data["M"].subs(v.P, p)
    n = data["N"]
    lift = data["lift"].subs(v.P, p)
    projector = data["complement_projector"].subs(v.P, p)
    assert projector.rank() == 56
    j = sp.Matrix([sp.Rational((7*i) % 11-5, 13) for i in range(60)])
    unconstrained = -m.inv()*j
    prescribed = sp.Matrix([2, -3, 5, 7])
    difference = prescribed-n*unconstrained
    constrained = unconstrained+lift*difference
    assert n*constrained == prescribed
    assert projector.T*(m*constrained+j) == sp.zeros(60, 1)
    action = (constrained.T*m*constrained)[0]/2+(j.T*constrained)[0]
    original = (unconstrained.T*m*unconstrained)[0]/2+(j.T*unconstrained)[0]
    assert action-original == (difference.T*data["D_inverse"].subs(v.P, p)*difference)[0]/2


def test_continuous_D_sign_proof_uses_actual_tube_not_too_wide_p_lower():
    b = v.domain_bounds()
    assert b["lower_p_squared_margin"] == sp.Rational(41, 10000) > 0
    assert b["upper_p_squared_margin"] == sp.Rational(11, 400) > 0
    assert b["polynomial_lower"] == sp.Rational(15623, 10**6) > 0
    assert sp.Poly(b["polynomial_derivative"], v.P).all_coeffs() == [3, 4, 1]
    assert b["D00_lower"] == sp.Rational(46869, 1100000) > 0
    assert b["minus_Dii_lower"] == sp.Rational(14, 55) > 0
    assert b["D_inverse_norm_upper"] == sp.Rational(1100000, 46869) < 24
    # The old quotient bound used .45 safely, but that bracket cannot prove THIS sign.
    assert b["polynomial"].subs(v.P, sp.Rational(9, 20)) < 0


def test_full_covariant_basis_and_CD_coefficients():
    assert _zero(v.rest_coefficients()["reconstruction_residual"])
    coefficients = v.cd_coefficients()
    assert all(value == 0 for value in coefficients["residuals"].values())
    assert coefficients["clock_values"] == {
        "a": 0, "b": -5/(2*v.H), "c": 0, "d": -1/v.H}
    assert coefficients["clock_X_derivatives"]["a"] == 15*v.HP/(8*v.H**2)
    assert coefficients["clock_X_derivatives"]["c"] == 1/v.H


def test_background_is_preserved_by_full_Euler_source_not_only_vector_projection():
    data = v.background()
    assert data["source"].shape == (64, 1)
    assert data["stationary"].shape == (64, 1)
    assert _zero(data["source"]) and _zero(data["stationary"])
    assert data["qx_clock"] == 0
    assert all(data[key] == 0 for key in ("c_clock", "quartic_clock", "p_phi_clock", "f_phi_clock"))
    assert sp.factor(sp.diff(data["a"], data["u"])/data["a"]-data["Hubble"]) == 0


def test_ADM_hessian_includes_shift_and_scale_before_the_vector_cancellation():
    data = v.adm_linearization()
    z = data["symbols"]
    eps, a, hubble = z["eps"], z["a"], z["Hubble"]
    h = data["Hessian"]
    assert h[0, 0] == -eps*z["n_dot"]
    for i in range(3):
        assert sp.expand(h[0, i+1]+eps*(z["n_gradient"][i]+hubble*z["shift"][i])) == 0
        # The g^{0i} term is required to remove Hubble*shift from H.v.
        assert sp.expand(data["Hv"][i+1]-eps*z["n_gradient"][i]) == 0
        for j in range(3):
            expected = (-a**2*hubble*int(i == j)
                        +eps*((z["shift_gradient"][i, j]+z["shift_gradient"][j, i])/2
                              -a**2*(2*hubble*z["zeta"]+z["zeta_dot"]-2*hubble*z["n"])*int(i == j)))
            assert sp.expand(h[i+1, j+1]-expected) == 0
    expected_box = (-3*hubble+eps*(z["n_dot"]-3*z["zeta_dot"]+6*hubble*z["n"]
                                  +sp.trace(z["shift_gradient"])/a**2))
    assert sp.expand(data["Box"]-expected_box) == 0


def test_coefficient_background_terms_are_not_discarded():
    data = v.adm_linearization()
    z = data["symbols"]
    independent = -3*z["n_dot"]/(2*v.H)+(15*v.HP/(4*v.H**2)-6*z["Hubble"]/v.H)*z["n"]
    assert sp.factor(data["first_variation"][0]-independent) == 0
    assert sp.factor(data["first_variation"][0]+3*z["n_dot"]/(2*v.H)) != 0
    assert _zero(data["first_variation_residual"])
    forbidden = {z["zeta"], z["zeta_dot"], *z["shift"], *z["shift_dot"],
                 *z["shift_gradient"], *z["zeta_gradient"]}
    assert not any(value.free_symbols & forbidden for value in data["first_variation_CD"])


def test_actual_center_first_variation_with_no_metric_jet_assumption():
    data = v.adm_linearization()
    z = data["symbols"]
    expected = sp.Matrix([-sp.Rational(3, 2)*z["n_dot"],
                          *(-sp.Rational(5, 2)*value for value in z["n_gradient"])])
    assert data["first_variation_CD"].subs({v.H: 1, z["Hubble"]: 0}) == expected
    assert v.schur()["D_inverse"].subs(v.P, sp.Rational(1, 2)) == sp.diag(8, -8, -8, -8)/3


def test_off_center_exact_gradient_shift_and_wrong_partial_shift_control():
    data = v.rolling()
    u, h, hubble = data["u"], data["h"], data["Hubble"]
    a = data["gradient_amplitude"]
    assert sp.factor(-sp.diff(a, u)-data["Vstar_lapse_value"]-21*hubble/(8*h)) == 0
    assert data["mass_spatial_shift"] == 1/h
    assert sp.factor(data["Vstar_curl_lapse_value"]-33*hubble/(8*h)) == 0
    # V+A(u)dn is not the claimed exact-gradient change when A' != 0.
    assert sp.diff(a, u).subs(u, sp.Rational(1, 2)) != 0
    position = sp.Symbol("position", real=True)
    lapse = sp.Function("n")(u, position)
    potential = a*lapse
    assert sp.diff(sp.diff(potential, position), u)-sp.diff(sp.diff(potential, u), position) == 0
    wrong = sp.diff(a*sp.diff(lapse, position), u)-sp.diff(a*sp.diff(lapse, u), position)
    assert sp.expand(wrong-sp.diff(a, u)*sp.diff(lapse, position)) == 0


def test_nonunit_physical_parameters_and_no_full_health_claim():
    physical = v.units(mass_squared=3, time_scale=2, kinetic_coefficient=5)
    assert physical["lambda_normalized"] == sp.Rational(5, 12)
    assert physical["isolated_center_mass_squared_physical"] == sp.Rational(8, 5)
    assert physical["isolated_center_mass_squared_normalized"] == sp.Rational(32, 5)
    assert physical["D_physical_factor"] == sp.Rational(1, 3)
    assert physical["V_normalized_factor"] == 2
    assert v.calibration()["full_coupled_health_gap_or_UV_claim"] is False
    assert v.calibration()["point_transformation_removes_lapse_velocity"] is False


@pytest.mark.parametrize("p", [sp.Rational(1, 2), sp.sqrt(sp.Rational(9, 40)), sp.sqrt(sp.Rational(11, 40))])
def test_closed_domain_accepts_exact_endpoints(p):
    assert v.require_domain(p) == (p, 1, 1)


@pytest.mark.parametrize("p", [True, 0.5, "1/2", 0, -1, sp.Rational(1, 4),
                                sp.Rational(3, 4), sp.sqrt(2)/4, sp.oo, sp.nan,
                                sp.I, sp.Symbol("p"), sp.Float("0.5"), [1]])
def test_domain_rejects_inexact_or_outside_chart(p):
    with pytest.raises((TypeError, ValueError)):
        v.require_domain(p)


@pytest.mark.parametrize("parameters", [{"mass_squared": 0}, {"mass_squared": True},
                                         {"time_scale": -1}, {"time_scale": "2"},
                                         {"kinetic_coefficient": 0}, {"kinetic_coefficient": 1.0},
                                         {"kinetic_coefficient": sp.I}])
def test_unit_guards(parameters):
    with pytest.raises((TypeError, ValueError)):
        v.units(**parameters)
