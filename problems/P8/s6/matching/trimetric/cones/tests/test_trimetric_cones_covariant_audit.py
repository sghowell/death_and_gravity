"""Independent coframe, physical-clock, source and full-solution checks."""

import pytest
import sympy as sp
from p8_trimetric import matter, model
from p8_trimetric_cones import background, bounds, rolling, tensor


def zero(value):
    values = list(value) if isinstance(value, sp.MatrixBase) else [value]
    assert all(sp.simplify(item) == 0 for item in values)


def test_all_coframe_components_give_the_actual_frame_null_cone_identity():
    d = background.derive()
    ne, ae, nv, av, nu, au = [d[key] for key in ("n_e", "a_e", "n_v", "a_v", "n_u", "a_u")]
    e, v, u = sp.diag(ne, ae, ae, ae), sp.diag(nv, av, av, av), sp.diag(nu, au, au, au)
    # Literal covariant matter gradient of a homogeneous isotropic source.
    ju = sp.diag(-au**3*d["rho"], *([nu*au**2*d["p"]]*3))
    eu = model.euler_maps(e, v, u, pg=d["p_g"], pf=d["p_f"], b=d["B"],
                          matter_gradient=ju, epsilon=d["epsilon"])["E_u"]
    lapse = eu[0, 0]/(-2*au**3)
    zero(lapse-d["lapse_residual"])
    for index in (1, 2, 3):
        space = eu[index, index]/(-2*nu*au**2)
        zero(space-d["spatial_residual"])
        zero(space-lapse-d["weighted_cone_residual"])


def test_weighted_cone_bound_is_a_positive_weight_average_not_two_separate_bounds():
    pe, pv, ce, cv = sp.symbols("P_e P_v c_e c_v", positive=True)
    average = (pe*ce+pv*cv)/(pe+pv)
    zero(average-1-(pe*(ce-1)+pv*(cv-1))/(pe+pv))
    # One mode can be subluminal even when the weighted excess is positive.
    assert average.subs({pe: 1, pv: 3, ce: sp.Rational(1, 2), cv: 2}) > 1


def test_nondiagonal_TT_matrix_second_variation_has_the_unit_polarization_factor():
    # An off-diagonal spatial tensor with norm²=28, unlike the diagonal
    # unit-polarization fixture in the primary implementation.
    shear = sp.zeros(4)
    shear[1:4, 1:4] = sp.Matrix([[0, 1, 2], [1, 0, 3], [2, 3, 0]])
    q, big_a, big_n = sp.Integer(3), sp.Integer(2), sp.Rational(2, 7)
    u = sp.diag(big_n, big_a, big_a, big_a)
    inv = u.inv()
    ge, gv, gu = sp.Integer(1), sp.Integer(-2), sp.Integer(3)
    up, upp = u*gu*shear/2, u*gu**2*shear**2/4

    def trace_second(amplitude):
        ep, epp = amplitude*shear/2, amplitude**2*shear**2/4
        return sp.trace(2*inv*up*inv*up*inv-inv*upp*inv-2*inv*up*inv*ep+inv*epp)

    # Exponential TT paths keep all three determinants fixed to this order;
    # homogeneous canonical matter has no time-gradient change along them.
    zero(sp.trace(inv*upp)+sp.trace(inv*up)**2-sp.trace((inv*up)**2))
    coefficient = -q*u.det()*(trace_second(ge)+trace_second(gv))
    expected = -u.det()*q/(4*big_a)*((ge-gu)**2+(gv-gu)**2)*sp.trace(shear**2)
    zero(coefficient-expected)
    assert coefficient == -696


def test_full_auxiliary_source_completion_keeps_local_contact_and_common_source():
    d = tensor.derive()
    pe, pv, ge, gv, gu, j = [d[key] for key in ("P_g", "P_f", "gamma_e", "gamma_v", "gamma_u", "j")]
    lag = -(pe*(ge-gu)**2+pv*(gv-gu)**2)/4+j*gu
    solution = sp.solve(sp.diff(lag, gu), gu)[0]
    zero(solution-d["sourced_auxiliary_map"])
    reduced = sp.cancel(lag.subs(gu, solution))
    zero(sp.diff(reduced, j)-solution)
    zero(sp.diff(reduced, j, 2)-2/(pe+pv))
    zero(sp.diff(reduced, j, ge)-pe/(pe+pv))
    zero(sp.diff(reduced, j, gv)-pv/(pe+pv))


def test_clock_volume_and_spatial_units_independently_recover_both_tensor_coefficients():
    d = tensor.derive()
    n, a, big_n, big_a, g = [d[key] for key in ("n", "a", "N", "A", "G")]
    measure_ratio = n*a**3/((big_n*n)*(big_a*a)**3)
    # Einstein metric proper derivatives have D_r=N D_h and k_r=A k_h.
    kinetic = 2*g*measure_ratio*big_n**2
    gradient = 2*g*measure_ratio*big_a**2
    zero(kinetic-d["G_T_h"])
    zero(gradient-d["F_T_h"])
    zero(gradient/kinetic-d["common_speed_squared"])


def test_unequal_Einstein_coefficients_do_not_inherit_exact_even_odd_decoupling():
    d = tensor.derive()
    cdot, rdot = sp.symbols("cdot rdot", real=True)
    density = (d["G"]*(cdot+rdot)**2+d["F"]*(cdot-rdot)**2)/8
    zero(sp.diff(density, cdot, rdot)-(d["G"]-d["F"])/4)
    assert sp.diff(density, cdot, rdot).subs({d["G"]: 2, d["F"]: 1}) != 0


def test_actual_rolling_family_solves_every_full_coframe_Euler_component():
    d = rolling.derive()
    a, n, q, g, eps = [d[key] for key in ("A", "N", "q", "G", "epsilon")]
    u = sp.diag(n, a, a, a)
    source = matter.canonical_source(u, [d["psi_prime"], 0, 0, 0])
    # B-sector R=-6(Hdot+2H²), hence G00=+3H². This is not A's FK sign.
    einstein = sp.diag(3*d["H_r_squared"], *([-2*d["H_r_prime"]-3*d["H_r_squared"]]*3))
    full = model.euler_maps(sp.eye(4), sp.eye(4), u, b=-6*q, pg=q, pf=q, bg=-q, bf=-q,
                            g=g, f=g, epsilon=eps, matter_gradient=source["J_u"],
                            einstein_g=einstein, einstein_f=einstein)
    for key in ("E_e", "E_v", "E_u"):
        zero(full[key])


def test_actual_scalar_current_and_physical_clock_do_not_use_a_fictitious_static_point():
    d = rolling.derive()
    x, a, n, h, adot, rho = [d[key] for key in ("x", "A", "N", "H_r", "A_prime", "rho")]
    current_log_rate = 3*(h+adot/a)+sp.diff(rho, x)*adot/(2*rho)
    zero(current_log_rate)
    zero((h+adot/a)/n-d["H_h"])
    positive_hubble = sp.sqrt(sp.expand(d["H_r_squared"]))/a
    zero(d["H_h"]-positive_hubble)
    assert positive_hubble.is_positive is True
    assert d["physical_speed"].subs(x, 1) == 7
    assert sp.limit(d["physical_speed"], x, 0, dir="+") == 1


def test_mixed_sign_auxiliary_point_is_not_mislabeled_as_a_full_flat_solution():
    e, v, u = sp.diag(sp.Rational(1, 2), 1, 1, 1), sp.eye(4), sp.eye(4)
    source = matter.canonical_source(u, [sp.sqrt(2), 0, 0, 0])
    result = model.euler_maps(e, v, u, pg=-2, pf=1, b=sp.Rational(5, 2),
                              matter_gradient=source["J_u"], epsilon=1)
    zero(result["E_u"])
    assert result["E_e"] != sp.zeros(4)
    assert result["E_v"] != sp.zeros(4)


def test_mixed_sign_full_vacuum_has_positive_spring_from_literal_eliminated_action():
    variation = sp.Symbol("variation", real=True)
    shear = sp.diag(0, 1, -1, 0)
    e, v = sp.eye(4)+variation*shear/2, sp.eye(4)-variation*shear/2
    u = 2*e-v  # -3*(pg*e+pf*v)/B at pg=-2,pf=1,B=3.
    lag = model.potential_density(e, v, u, b=3, pg=-2, pf=1, bg=2, bf=-1)
    coefficient = sp.diff(lag, variation, 2).subs(variation, 0)/2
    # tr(shear²)=2, gamma_e-gamma_v=2*variation, density=-q_eff*variation².
    assert coefficient == -4
    assert background.controls()["mixed_sign_full_flat_relative_spring"] == -coefficient


def test_remainder_budget_uses_actual_normalized_null_stress_not_epsilon_alone():
    a, q, eps, bare = sp.symbols("A q epsilon bare", positive=True)
    actual = eps*bare
    data = bounds.speed_budget(a, q, actual, 1)
    fixed_actual = sp.Symbol("fixed_actual", positive=True)
    zero(data["speed_gap"].subs(bare, fixed_actual/eps)-a*fixed_actual/(4*q))
    normalized = bounds.speed_budget(2, 1, 12, sp.Rational(1, 14), sp.Rational(1, 56), sp.Rational(1, 56))
    assert normalized["corrected_kinetic_lower_bound"] > 0
    assert normalized["corrected_F_minus_G_lower_bound"] == sp.Rational(95, 28)


def test_asymmetric_high_frequency_cone_does_not_decide_a_light_only_effective_cone():
    # A positive quadratic control only, not a complete trimetric solution.
    momentum_squared, speed_squared, spring = sp.symbols("k2 speed2 spring", positive=True)
    kinetic = sp.diag(100, 1)
    gradient = sp.diag(25, 4)
    mass = spring*sp.Matrix([[1, -1], [-1, 1]])
    determinant = (momentum_squared*speed_squared*kinetic-momentum_squared*gradient-mass).det()
    low_energy = sp.solve(sp.diff(determinant, momentum_squared).subs(momentum_squared, 0), speed_squared)
    assert low_energy == [sp.Rational(29, 101)]
    assert (kinetic.inv()*gradient).eigenvals() == {sp.Rational(1, 4): 1, sp.Integer(4): 1}


@pytest.mark.parametrize("invalid", [0, -1, 0.5, True])
def test_cone_budget_does_not_admit_singular_wrong_sign_or_inexact_link(invalid):
    with pytest.raises((TypeError, ValueError)):
        bounds.speed_budget(1, invalid, 1, 1)
