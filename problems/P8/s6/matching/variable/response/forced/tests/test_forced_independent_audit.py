"""Root audit: literal physical action, moving maps and causal norm proof.

The reference physical equations below are reconstructed from the two
Einstein TT kinetic terms and the actual f-null source, not the candidate's
canonical operator. Exact fixtures audit moving-map jets; the written
proof and coefficientwise bounds, not these fixtures, certify the domain.
"""

from fractions import Fraction
from functools import cache

import pytest
import sympy as sp
from p8_variable_forced import bounds, coefficients, kernel, operator, source

F = Fraction


@cache
def literal_action():
    u = sp.Symbol("audit_u", real=True)
    delta = sp.Symbol("audit_delta", positive=True)
    momentum = sp.Symbol("audit_K", nonnegative=True)
    c, d = 2+delta, 1+u**2
    a, b = d**2, 2/d**2
    y, h = b/a, sp.diff(a, u)/a
    kg, kf = a**3/8, b**3/(8*c)
    gg, gf = a/8, c*b/8
    # -2 D_f H_f=(c-y)P/(c*y^3), H_f=-h/c for constant c.
    link = 2*y**3*sp.diff(h, u)/(c*(c-y))
    spring = a**3*y*link/8
    mass_matrix = sp.diag(kg, kf)
    potential = sp.Matrix([[gg*momentum+spring, -spring],
                           [-spring, gf*momentum+spring]])
    drag = mass_matrix.inv()*mass_matrix.diff(u)
    stiffness = mass_matrix.inv()*potential
    total = kg+kf
    fs = sp.sqrt(2*total)
    fr = sp.sqrt(2*kg*kf/total)
    position = sp.Matrix([[1/fs, -(kf/total)/fr],
                          [1/fs, (kg/total)/fr]]).applyfunc(sp.simplify)
    # Variation of +(a^3/2)*sigma*gamma_g, including the factor 2 in
    # variation of Kg*gamma_g'^2, gives the literal force (2*sigma,0).
    forcing = sp.Matrix([a**3/(4*kg), 0])
    return {"u": u, "delta": delta, "K": momentum, "a": a,
            "kinetic": mass_matrix, "drag": drag, "stiffness": stiffness,
            "M": position, "Mprime": position.diff(u),
            "Msecond": position.diff(u, 2), "force": forcing}


FIXTURES = [(sp.S.Zero, sp.Rational(1, 10**9), 0),
            (sp.Rational(1, 100), sp.Rational(1, 10**12), 4),
            (-sp.Rational(1, 200), sp.Rational(1, 10**10), 1)]


def zero_matrix(matrix):
    assert all(sp.simplify(value) == 0 for value in matrix)


@pytest.mark.parametrize("time,delta,momentum", FIXTURES)
def test_literal_action_with_both_moving_normalization_jets(time, delta, momentum):
    ref, actual = literal_action(), operator.derive()
    rv = {ref["u"]: time, ref["delta"]: delta, ref["K"]: momentum}
    av = {actual["u"]: time, actual["delta"]: delta, actual["K"]: momentum}
    at = lambda matrix: matrix.subs(rv).applyfunc(sp.simplify)
    mm, mp, mpp = (at(ref[key]) for key in ("M", "Mprime", "Msecond"))
    kinetic, drag, stiffness = (at(ref[key]) for key in ("kinetic", "drag", "stiffness"))
    zero_matrix(mm.T*kinetic*mm-sp.eye(2)/2)
    first = mm.inv()*(2*mp+drag*mm)
    zeroth = mm.inv()*(mpp+drag*mp+stiffness*mm)
    target_first = sp.Matrix([[0, actual["u"]*actual["dc"]],
                              [actual["u"]*actual["fc"], 0]]).subs(av)
    target_zeroth = sp.Matrix([[actual["A"], actual["C"]],
                               [actual["E"], actual["B"]]]).subs(av)
    target_force = sp.Matrix([actual["jL"], actual["jH"]]).subs(av)
    zero_matrix(first-target_first)
    zero_matrix(zeroth-target_zeroth)
    zero_matrix(mm.inv()*at(ref["force"])-target_force)
    zero_matrix(mm[0, :]-sp.Matrix([[actual["ag"], actual["bg"]]]).subs(av))


def test_center_second_jet_cannot_be_frozen_away():
    ref = literal_action()
    values = {ref["u"]: 0, ref["delta"]: sp.Rational(1, 10**9), ref["K"]: 0}
    zero_matrix(ref["Mprime"].subs(values))
    assert any(sp.simplify(value) != 0 for value in ref["Msecond"].subs(values))
    assert ref["force"] == sp.Matrix([2, 0])


@pytest.mark.parametrize("polarization", [sp.diag(1, -1, 0),
                                         sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]])])
def test_external_spatial_TT_probe_is_background_covariantly_conserved(polarization):
    time, z, momentum = sp.symbols("t z k", real=True)
    scale, stress = sp.Function("a")(time), sp.Function("Pi")(time)
    contravariant = stress*sp.cos(momentum*z)*polarization/scale**2
    # For ds^2=dt^2-a^2 dx^2 and T^{0mu}=0, the only connection
    # contribution is Gamma^0_ij T^{ij}=a*a' trace(T^{ij}).
    temporal = scale*sp.diff(scale, time)*sp.trace(contravariant)
    spatial = contravariant[2, :].diff(z)
    assert temporal == 0
    zero_matrix(spatial)
    zero_matrix(spatial.subs(momentum, 0))
    assert polarization == polarization.T


def test_off_diagonal_derivative_is_the_formal_adjoint_not_a_symmetric_potential():
    d = operator.derive()
    u = d["u"]
    ell, heavy = sp.Function("ell")(u), sp.Function("heavy")(u)
    b_ell = d["E"]*ell+u*d["fc"]*sp.diff(ell, u)
    adj_h = (d["E"]-sp.diff(u*d["fc"], u))*heavy-u*d["fc"]*sp.diff(heavy, u)
    assert sp.simplify(heavy*b_ell-ell*adj_h-sp.diff(u*d["fc"]*ell*heavy, u)) == 0
    assert sp.simplify(d["C"]-d["E"]+sp.diff(u*d["fc"], u)) == 0
    assert sp.simplify((d["C"]-d["E"]).subs(u, 0)) != 0


def test_Liouville_transform_from_general_clock_jets():
    u, eps = sp.symbols("u eps", real=True, positive=True)
    rho = sp.sqrt(u*u+eps)
    amp = sp.sqrt(rho)
    zprime = 1/rho
    # Differentiate amp*psi(z(u)) coefficient by coefficient before
    # any substitution, independently of the candidate's chain rule.
    coeff_second = sp.simplify(rho**sp.Rational(3, 2)*amp*zprime**2)
    coeff_first = sp.simplify(rho**sp.Rational(3, 2)*(2*sp.diff(amp, u)*zprime+amp*sp.diff(zprime, u)))
    coeff_zero = sp.simplify(rho**sp.Rational(3, 2)*(sp.diff(amp, u, 2)+10*amp/rho**2))
    assert coeff_second == 1
    assert coeff_first == 0
    assert sp.simplify(coeff_zero-sp.Rational(39, 4)-3*eps/(4*rho**2)) == 0
    assert sp.simplify(rho*zprime-1) == 0


def test_single_Jost_solution_and_positive_retarded_jump_algebraically():
    x, mu = sp.symbols("x mu", real=True, positive=True)
    ff, fx = sp.symbols("F Fx")
    # x=(1-tanh(z))/2, x_z=-2x(1-x). Start from the standard
    # Gauss equation x(1-x)Fxx+(1-i*mu-2x)Fx+3F/4=0.
    xz = -2*x*(1-x)
    xzz = sp.diff(xz, x)*xz
    fxx = (-(1-sp.I*mu-2*x)*fx-sp.Rational(3, 4)*ff)/(x*(1-x))
    residual = xz**2*fxx+(xzz+2*sp.I*mu*xz)*fx+3*x*(1-x)*ff
    assert sp.expand(residual) == 0
    assert sp.simplify((sp.I*mu-(-sp.I*mu))/(2*sp.I*mu)) == 1
    assert kernel.MU**2 == sp.Rational(39, 4)
    # This infinity normalization only fixes an auxiliary ODE basis.
    # It does not prescribe an asymptotic state of the physical system.


def test_continuous_pole_bound_without_sampling_or_delta_derivatives():
    v, delta = sp.symbols("v delta", nonnegative=True)
    denominator = (2+delta)*(1+v)**4-2
    base = delta+8*v
    remainder = sp.Poly(sp.expand(denominator-base), delta, v)
    assert all(value > 0 for value in remainder.coeffs())
    assert sp.expand(base**2-32*delta*v-(delta-8*v)**2) == 0
    # Integrate the separate derivative bounds along a rectangle from
    # (delta,v)=(0,0). D>=delta+8v controls both contributions.
    assert max(F(32), F(1008, 8)) == 126
    dv_bound = 4+F(6, 100)+F(4, 10000)+F(1, 10**6)
    v2_bound = 12+F(8, 100)+F(2, 10000)
    assert dv_bound < F(41, 10) and v2_bound < F(121, 10)
    total = 126+80*(F(41, 320)+F(121, 640))
    assert total == F(1211, 8) < 152
    assert total+9 < 161
    assert coefficients.envelope()["mass_remainder"] == sp.Rational(total)


@pytest.mark.parametrize("radius", [F(1, 100), F(1, 1000), F(1, 10**6)])
def test_rebuild_causal_error_chain_with_Fraction_operator_norms(radius):
    r, ratio, mu_lower = radius, F(101, 100), F(3)
    rho_upper = ratio*r
    eta = F(161)*2*r*rho_upper/mu_lower
    exponent = F(1, 2)+eta
    # n!>=2*3^(n-2), n>=2, majorizes exp(exponent) by this rational.
    assert 1+exponent+exponent**2/(2*(1-exponent/3)) < 2
    heavy = F(4, 3)*r*rho_upper
    weighted_derivative = F(14, 3)*r*rho_upper
    adjoint = 12*heavy+10*weighted_derivative
    assert adjoint < 64*r*r
    # Independent resolvent identity GH-G0=-G0*b_delta*GH,
    # using the same proved norm bound for both retarded inverses.
    kernel_difference = heavy*161*heavy
    assert kernel_difference < 292*r**4
    light = 4*r*r/(1-18*r*r)
    assert light < 5*r*r
    bop = 10+60*r*r
    feedback = 5*r*r*64*r*r*bop
    assert feedback < F(1, 10000)
    full_light = 5*r*r*(1+64*r*r)/(1-feedback)
    assert full_light < 6*r*r
    light_error = 5*r*r*64*r*r*(1+bop*6*r*r)
    assert light_error < 323*r**4
    approximate_force = 1+bop*5*r*r
    assert approximate_force < F(101, 100)
    heavy_error = F(292)*r**4*F(101, 100)+heavy*bop*323*r**4
    assert heavy_error < 296*r**4
    assert (323+2*296)*r**4 < 1000*r**4
    reported = bounds.calibration(sp.Rational(r))
    assert reported["GH_over_r2"] == sp.Rational(heavy/(r*r))
    assert reported["GH_minus_G0_over_r4"] == sp.Rational(kernel_difference/r**4)
    assert reported["feedback"] == sp.Rational(feedback)
    assert reported["l_error_over_r4"] == sp.Rational(light_error/r**4)


def test_projection_source_is_derived_from_forced_symplectic_coefficients():
    ge, go, gep, gop, aa, omega, sigma = sp.symbols("ge go gep gop a Omega sigma", real=True)
    alpha_prime = -aa**3*go*sigma/(2*omega)
    beta_prime = aa**3*ge*sigma/(2*omega)
    assert sp.expand(alpha_prime*ge+beta_prime*go) == 0
    slope = sp.expand(alpha_prime*gep+beta_prime*gop)
    wronskian = ge*gop-gep*go
    assert sp.simplify(slope-source.projected_source_factor(aa, wronskian, omega)*sigma) == 0
    assert F(1)*F(4, 5)/(2*F(1)) == F(2, 5)
    assert F(2)-F(2, 5) == F(8, 5)


def test_compact_pulse_RMS_diagnostic_has_the_correct_direction():
    r = F(1, 100)
    # A lower stiffness estimate proves positivity; an UPPER endpoint
    # estimate bounds inf sqrt(B) above for the frequency/proxy ratio.
    lower = F(80)/(8+F(1, 100))-161*r*r
    endpoint_upper = 10+161*r*r
    assert lower > 9
    assert endpoint_upper < F(13, 4)**2
    assert F(3, 2)/F(13, 4) == F(6, 13)
    actual = bounds.calibration(sp.Rational(r))
    assert actual["stiffness_lower_times_r2"] == sp.Rational(lower)
    assert actual["stiffness_endpoint_upper_times_r2"] == sp.Rational(endpoint_upper)


def test_retarded_support_and_no_false_fixed_radius_delta_convergence():
    delta = sp.Rational(1, 10**12)
    assert kernel.universal_kernel(-sp.Rational(1, 200), 0, delta, momentum_squared=0) == 0
    for momentum in (0, 1, 4):
        assert bounds.physical_response_error(delta, 3, momentum_squared=momentum) == sp.Rational(3, 10**5)
        assert bounds.physical_response_error(delta**2, 3, momentum_squared=momentum) == sp.Rational(3, 10**5)
    with pytest.raises(ValueError):
        bounds.physical_response_error(delta, 1, sp.Rational(1, 10**6))
