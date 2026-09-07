"""Separately authored coframe/source and geometric audit.

The source is differentiated from its elementary-polynomial generating
function, not from the primary stress formulas. The control ODEs and
center points below are explicitly not full reconstructed parent solutions.
"""

from fractions import Fraction as Q
from functools import cache
from math import comb

import pytest
import sympy as sp
from p8_reciprocal_geometry import identities
from p8_variable_beta import model


@cache
def _literal_source():
    a, b, ng, nf, t = sp.symbols("a b Ng Nf t", positive=True)
    beta = sp.symbols("beta0:5", real=True)
    generating = sp.Poly((1+t*nf/ng)*(1+t*b/a)**3, t)
    density = sp.expand(-2*ng*a**3*sum(beta[j]*generating.nth(j) for j in range(5)))
    rho = -sp.diff(density, nf)/b**3
    pressure = sp.diff(density, b)/(3*nf*b**2)
    return a, b, ng, nf, beta, density, rho, pressure


def test_full_f_coframe_variations_from_all_five_elementary_coefficients():
    a, b, ng, nf, beta, density, rho, pressure = _literal_source()
    y, c = b/a, nf/ng
    P = 2*(beta[1]+2*beta[2]*y+beta[3]*y*y)
    assert sp.factor(rho+pressure-(c-y)*P/(c*y**3)) == 0
    # beta0 is absent from both f variations; beta4 is a null-free cosmological term.
    assert sp.diff(rho, beta[0]) == sp.diff(pressure, beta[0]) == 0
    assert sp.diff(rho, beta[4]) == 2
    assert sp.diff(pressure, beta[4]) == -2
    assert sp.diff(density, nf, 2) == 0


def test_literal_source_matches_the_pinned_model_without_constant_clock_assumption():
    a, b, ng, nf, beta, _density, rho, pressure = _literal_source()
    d = model.derive()
    substitutions = {a: 1, b: d["y"], ng: 1, nf: d["c"]}
    substitutions.update(dict(zip(beta, d["betas"])))
    assert sp.factor(rho.subs(substitutions)-d["rho_f"]) == 0
    assert sp.factor(pressure.subs(substitutions)-d["pressure_f"]) == 0
    # Variation holds phi fixed, not beta'(phi)=0 along the solution.
    phi = sp.Symbol("phi", real=True)
    clock_coefficients = {beta[j]: sp.Function(f"b{j}")(phi) for j in range(5)}
    literal_clock = _density.subs(clock_coefficients)
    assert sp.factor(-sp.diff(literal_clock, nf)/b**3-rho.subs(clock_coefficients)) == 0
    assert sp.diff(literal_clock, phi) != 0


def _fraction_density(a, b, ng, nf, beta):
    # Direct expansion of Ng*a^3 product(1+t*root_eigenvalue), before variation.
    result = Q(0)
    for j in range(5):
        if j <= 3:
            result += beta[j]*comb(3, j)*ng*a**(3-j)*b**j
        if j >= 1:
            result += beta[j]*comb(3, j-1)*nf*a**(4-j)*b**(j-1)
    return -2*result


@pytest.mark.parametrize("a,b,ng,nf,beta", [
    (Q(3, 2), Q(5, 4), Q(2, 3), Q(7, 6), (Q(1, 3), Q(-2), Q(5, 7), Q(4, 3), Q(-9, 5))),
    (Q(1), Q(2), Q(1), Q(4), (Q(-26), Q(4), Q(0), Q(0), Q(-1, 2))),
    (Q(7, 3), Q(4, 5), Q(5, 6), Q(9, 4), (Q(2), Q(3), Q(-4), Q(5), Q(7))),
])
def test_independent_fraction_polynomial_variation(a, b, ng, nf, beta):
    step_n, step_b = nf/10, b/10
    lapse_derivative = (_fraction_density(a, b, ng, nf+step_n, beta)
                        -_fraction_density(a, b, ng, nf-step_n, beta))/(2*step_n)
    # The exact four-node derivative stencil differentiates this cubic in b.
    scale_derivative = (_fraction_density(a, b-2*step_b, ng, nf, beta)
                        -8*_fraction_density(a, b-step_b, ng, nf, beta)
                        +8*_fraction_density(a, b+step_b, ng, nf, beta)
                        -_fraction_density(a, b+2*step_b, ng, nf, beta))/(12*step_b)
    rho, pressure = -lapse_derivative/b**3, scale_derivative/(3*nf*b**2)
    y, c = b/a, nf/ng
    P = 2*(beta[1]+2*beta[2]*y+beta[3]*y*y)
    assert rho+pressure == (c-y)*P/(c*y**3)


def test_lapse_first_einstein_scale_equation_retains_general_positive_F():
    _a, b, _ng, nf, _beta, potential, rho, pressure = _literal_source()
    F, bd, bdd, nfd = sp.symbols("F bdot bddot Nfdot", positive=True)
    lagrangian = -3*F*b*bd**2/nf+potential
    dt_velocity_derivative = (bd*sp.diff(sp.diff(lagrangian, bd), b)
                              +bdd*sp.diff(lagrangian, bd, 2)
                              +nfd*sp.diff(sp.diff(lagrangian, bd), nf))
    lapse = sp.diff(lagrangian, nf)/b**3
    scale = (sp.diff(lagrangian, b)-dt_velocity_derivative)/(3*nf*b**2)
    Hf = bd/(nf*b)
    DfHf = bdd/(nf**2*b)-bd*nfd/(nf**3*b)-bd**2/(nf**2*b**2)
    assert sp.factor(lapse-(3*F*Hf**2-rho)) == 0
    assert sp.factor(scale-(2*F*DfHf+3*F*Hf**2+pressure)) == 0
    assert sp.factor(scale-lapse-(2*F*DfHf+rho+pressure)) == 0


def test_unfixed_lapse_reciprocal_raychaudhuri_and_off_shell_integrating_factor():
    u = sp.Symbol("u", real=True)
    a, c = sp.Function("a", positive=True)(u), sp.Function("c", positive=True)(u)
    alpha, F, P = sp.symbols("alpha F P", positive=True)
    b, y = alpha/a, alpha/a**2
    h = sp.diff(a, u)/a
    z, Z = h/c, h/y
    D = Z-z
    Hf = sp.diff(b, u)/(c*b)
    assert sp.factor(Hf+z) == 0
    source = 2*F*sp.diff(z, u)-(c-y)*P/y**3
    equation = sp.diff(D, u)+P*D/(2*F*y**2*z)-sp.diff(Z, u)
    assert sp.factor(equation+source/(2*F)) == 0
    assert sp.diff(source, sp.diff(c, u)) != 0  # the lapse derivative was not frozen
    assert sp.factor(sp.diff(Z, u)-sp.diff(a*sp.diff(a, u), u)/alpha) == 0


def test_center_sign_does_not_divide_by_h_or_z():
    hp, c, y, F = sp.symbols("hp c y F", positive=True)
    P = 2*F*y**3*hp/(c*(c-y))
    zp, Dp = hp/c, hp/y-hp/c
    assert sp.factor(2*F*y**3*zp-(c-y)*P) == 0
    assert sp.factor(Dp-hp*(c-y)/(c*y)) == 0
    assert P.subs({hp: 4, c: 4, y: 2, F: 3}) == 24
    assert P.subs({hp: 4, c: 1, y: 2, F: 3}) == -192
    # This is an algebraic point control, not an actual sourced solution.
    f_null = 2*F*hp/c**2
    assert sp.factor(2*F*zp-c*f_null) == 0
    assert f_null.is_positive  # added f matter permits the c=y center outside this theorem


def test_integrating_factor_identity_and_geometric_hypothesis_control():
    u = sp.Symbol("u", real=True)
    D, integrating = sp.Function("D")(u), sp.Function("I", positive=True)(u)
    rate, forcing = sp.symbols("rate forcing", nonnegative=True)
    residual = sp.diff(integrating*D, u).subs(
        {sp.diff(integrating, u): rate*integrating, sp.diff(D, u): forcing-rate*D})
    assert sp.factor(residual-integrating*forcing) == 0
    # Removing B>=0 invalidates positivity even in the scalar comparison ODE.
    counter = sp.exp(-u)*(1-2*u)
    assert sp.factor(sp.diff(counter, u)+counter) == -2*sp.exp(-u)
    assert counter.subs(u, 0) == 1 and counter.subs(u, 1) < 0


def test_affine_primitive_is_a_finite_window_bound_with_fixed_initial_z():
    a, alpha, h, z0, dz = sp.symbols("a alpha h z0 dz", positive=True)
    z = z0+dz
    measure = alpha*h/(a*z)
    primitive_derivative = -alpha*h/(a*z0)
    gap = sp.factor(measure+primitive_derivative)
    assert gap == -alpha*h*dz/(a*z0*(z0+dz))
    assert gap.is_negative
    a0, a1 = Q(5, 4), Q(7, 2)
    alpha0, z00 = Q(2), Q(3, 2)
    finite_window = alpha0/z00*(1/a0-1/a1)
    assert finite_window == Q(24, 35)
    assert 0 < finite_window < alpha0/(z00*a0)
    primary = identities.derive()
    substitutions = dict(zip(primary["variables"], (a, alpha, h, 1, z, z0, 1)))
    assert primary["affine_tail_bound"].subs(substitutions) == alpha/(z0*a)


def test_geodesic_pullback_and_comoving_vs_noncomoving_scope():
    b, c, momentum, a0, growth = sp.symbols("b c momentum a0 growth", positive=True)
    null_du, dx = momentum/(b*c), momentum/b**2
    assert sp.factor(c**2*null_du**2-b**2*dx**2) == 0
    proper = c/sp.sqrt(1+momentum**2/b**2)
    assert sp.factor((proper/(b*c/momentum))**2) == momentum**2/(b**2+momentum**2)
    assert proper.subs(momentum, 0) == c
    assert sp.limit(proper/(b*c/momentum), momentum, 0, dir="+") == 0
    # g proper-time integrand stays uniformly positive for a>=a0; null one grows.
    g_integrand_squared = (a0+growth)**2/((a0+growth)**2+momentum**2)
    baseline = a0**2/(a0**2+momentum**2)
    assert sp.factor(g_integrand_squared-baseline).is_positive


def test_nonconstant_affine_primitive_omission_and_shift_vs_TT_scope():
    a, alpha, h, z, zp = sp.symbols("a alpha h z zp", positive=True)
    # A moving z is not a constant of integration; its derivative cannot be omitted.
    primitive = alpha/(a*z)
    dt = a*h*sp.diff(primitive, a)+zp*sp.diff(primitive, z)
    assert sp.factor(alpha*h/(a*z)+dt) == -alpha*zp/(a*z**2)
    y, c, beta1, beta2, beta3 = sp.symbols("y c beta1 beta2 beta3", real=True)
    shift_P = 2*(beta1+2*beta2*y+beta3*y**2)
    tensor_mu = 2*y*(beta1+beta2*(c+y)+beta3*c*y)
    assert sp.factor(tensor_mu-y*shift_P-2*y*(c-y)*(beta2+beta3*y)) == 0
    assert (tensor_mu-y*shift_P).subs({y: 2, c: 4, beta1: 0, beta2: 1, beta3: 0}) == 8
