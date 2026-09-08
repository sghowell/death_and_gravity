"""Exact nonlinear chart, physical scalar function and boundary controls."""
import sympy as sp
from p8_affine_aligned import alignment
from p8_affine_nonlinear import adm, constraints


def test_full_time_and_spatial_ADM_identities():
    data = adm.transform()
    for name in ("trace_velocity_residual", "lapse_velocity_squared_residual", "remaining_linear_V",
                 "full_time_action_residual", "full_spatial_lapse_gradient_residual"):
        assert data[name] == 0
    assert all(value == 0 for value in constraints.boundaries().values())
    assert not data["Lv"].has(adm.K, adm.KK, adm.V)


def test_complete_original_scalar_potential_and_lapse_Jacobian():
    data = adm.clock_lapse()
    assert data["actual_background_lapse_equation"] == 0
    assert data["actual_secondary_lapse_Jacobian"] == 0
    assert data["I_N"] == 0
    assert data["I_NN"] == -18*adm.u/(1+adm.u**2)**7


def test_omitting_linear_velocity_primitive_changes_the_actual_Jacobian():
    data, bg = adm.transform(), alignment.parent.old.background()
    N = adm.N
    volume, B, linear, potential = [data[name].subs(adm.s, 1/N)
                                    for name in ("volume", "B", "Blinear", "Fnew")]
    wrong = N*((-2*bg["H"]-volume*linear)**2/(4*sp.Rational(2, 3)*volume*B)
               -volume*potential+bg["ell"]**2/(2*volume))
    bad_Jacobian = sp.factor(sp.diff(wrong, N, 2).subs(N, 1))
    assert sp.factor(bad_Jacobian-adm.clock_lapse()["second"]) != 0
    assert (bad_Jacobian-adm.clock_lapse()["second"]).subs(adm.u, 0) == 18


def test_omitting_spatial_curvature_boundary_leaves_a_lapse_gradient():
    data = adm.coefficients()
    B4, delta = -data["f"], -2*adm.X*data["omega_X"]
    wrong = sp.factor(data["E"]-2*B4*delta**2)
    assert wrong != 0
    assert wrong.subs({adm.X: 1, adm.u: 0}) != 0
