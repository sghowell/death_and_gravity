"""Root-owned scaled full-ODE, limiting dynamics and physical-amplitude audit."""

import sympy as sp
from p8_composite_asymptotics import limiting, scaled


def test_direct_full_positive_root_ode_has_the_claimed_regular_scaling():
    data = scaled.core()
    e, eta, j = data["e"], data["eta"], data["j"]
    y, density, h = 1/e, eta*e**2, j*e
    xx, rr = data["X_bar"]/e, -e*data["R_bar"]
    yd = y*(1+y)*(xx*rr-h*(xx+rr))/(xx-y*rr)
    rd = -3*h*(density+2*y/(1+y)**2)
    ratio = (xx-(1+y)*h)/((1+y)*h-y*rr)
    assert sp.cancel(yd/y-data["lambda"]) == 0
    assert sp.cancel(-yd/y**2-data["eprime"]) == 0
    assert sp.cancel(rd/e**2-2*eta*data["eprime"]/e-data["etaprime"]) == 0
    assert sp.cancel(e*ratio-data["kappa"]) == 0
    assert data["null_over_e"].subs(e, 0) == 2


def test_regular_kinetic_templates_come_from_literal_physical_action():
    data = scaled.canonical_jets()
    e, k = data["e"], data["kappa"]
    ae, mass = sp.symbols("Ae M", positive=True)
    a_g, y, c = ae*e/(1+e), 1/e, k/e
    k1 = mass**2*a_g**3*(1+c)/8
    k2 = mass**2*a_g**3*y**3*(1+c)/(8*c)
    total, relative = k1+k2, k1*k2/(k1+k2)
    assert sp.cancel(8*total/(mass**2*ae**3)-data["shape_sum"]) == 0
    assert sp.cancel(8*relative/(mass**2*ae**3*e**2)-data["shape_relative"]) == 0
    mu = y*(y-1)*(y-c)/(1+y)
    algebraic = mu*(1+c/y**3)/(1+c)**2
    assert sp.cancel(algebraic-data["mass_squared"]) == 0


def test_exact_cd_family_clock_is_not_a_constant_acceleration_substitute():
    data = scaled.cd_family()
    u, eps, a0 = data["u"], data["epsilon"], data["a0"]
    target = 4*u/(4/(a0*eps)+u**2)
    assert sp.cancel(data["h"]-target) == 0
    assert sp.cancel(sp.diff(data["A_e"], u)/data["A_e"]-target) == 0
    assert sp.cancel(sp.diff(target, u)/eps).subs(eps, 0) == a0
    assert sp.cancel(sp.diff(target, u, 2)/eps).subs(eps, 0) == 0
    assert sp.diff(target, u, 2) != 0  # only the fixed-u leading limit is zero


def test_limiting_background_is_derived_from_both_scaled_equations():
    data = scaled.core()
    rr, j, ac = sp.symbols("R j ac", positive=True)
    eta = 3*(rr**2-1)
    at = {data["e"]: 0, data["eta"]: eta, data["j"]: j, data["a_current"]: ac}
    ep_over_e = -data["lambda"].subs(at)
    rprime = sp.simplify(data["etaprime"].subs(at)/(6*rr))
    jprime = sp.simplify(data["jprime"].subs(at))
    assert sp.simplify(ep_over_e-rr-j) == 0
    assert sp.simplify(rprime-1+rr*(rr+j)) == 0
    assert sp.simplify(jprime-ac+j*(rr+j)) == 0
    assert sp.simplify(rprime+jprime-1-ac+(rr+j)**2) == 0


def test_exact_forward_solution_and_bounds_do_not_use_finite_time_sampling():
    data = limiting.solution()
    u, z = data["u"], data["z"]
    assert sp.simplify(sp.diff(z, u, 2)-z-11) == 0
    assert z.subs(u, 0) == 1 and sp.diff(z, u).subs(u, 0) == 2
    s = sp.Symbol("exp_u", positive=True)
    z_s, zp_s = 7*s+5/s-11, 7*s-5/s
    # s>=1 is the exact forward-time domain, not a numerical sample.
    assert sp.factor(z_s-1-(s-1)*(7*s-5)/s) == 0
    assert sp.factor(zp_s-z_s-(11-10/s)) == 0
    assert sp.factor(zp_s**2-z_s**2-22*z_s+19) == 0
    zz = sp.Symbol("z", positive=True)
    speed2 = 1+22/zz-19/zz**2
    assert sp.factor(sp.Rational(140, 19)-speed2-(11*zz-19)**2/(19*zz**2)) == 0
    assert sp.Rational(140, 19) < 9
    assert sp.simplify(sp.diff(data["positive_R_numerator"], u)-z) == 0
    assert data["positive_R_numerator"].subs(u, 0) == 2
    assert z_s.subs(s, sp.Rational(9, 10)) < 1  # backward extrapolation fails


def test_full_time_relative_equation_has_a_positive_comparison_remainder():
    data = limiting.general_jets()
    v, a = data["v"], data["a"]
    tr = (3*v**2-1-a)/(2*v)
    derivative = lambda f: sp.diff(f, v)*(1+a-v**2)-sp.diff(f, a)*a*v
    curvature = sp.factor(derivative(tr)+tr**2)
    growth = curvature-(v**2-v)
    # Positive on the proved forward orbit: 1<=v<3 and a>0.
    margin = 3*(a+1)**2/(4*v**2)+(v-1)*(3-v)/4
    assert sp.factor(growth-sp.Rational(1, 4)-margin) == 0
    assert sp.factor(curvature-data["N_relative"]) == 0
    assert sp.factor(growth-data["heavy_growth_coefficient"]) == 0
    # The comparison is for the evolving coefficient, not its initial value.
    assert derivative(growth).subs({v: 2, a: 11}) != 0


def test_positive_green_comparison_identity_is_an_evolution_statement():
    t, s = sp.symbols("t s", nonnegative=True)
    h = sp.Function("H")
    coefficient = sp.Function("C")
    kernel = 2*sp.sinh((t-s)/2)
    integral = sp.Integral(kernel*(coefficient(s)-sp.Rational(1, 4))*h(s), (s, 0, t))
    trial = sp.cosh(t/2)+integral
    assert sp.simplify(sp.diff(trial, t, 2)-trial/4
                       -(coefficient(t)-sp.Rational(1, 4))*h(t)) == 0
    assert trial.subs(t, 0).doit() == 1
    assert sp.diff(trial, t).subs(t, 0).doit() == 0


def test_physical_composite_projection_and_small_amplitude_limit_are_distinct():
    data = scaled.canonical_jets()
    e, k = data["e"], data["kappa"]
    ll, hh = sp.symbols("l H", real=True)
    fs, fr = sp.symbols("f_sum f_relative", positive=True)
    w1 = k*e**2/(1+k*e**2)
    hg, hf = ll/fs-(1-w1)*hh/fr, ll/fs+w1*hh/fr
    physical = (e*hg+hf)/(1+e)
    normalized = physical.subs(fr, fs*e*sp.sqrt(k)/(1+k*e**2))
    assert sp.cancel(fs*normalized-ll-data["physical_source_relative"]*hh) == 0
    eps, z, v, mass = sp.symbols("epsilon z v M", positive=True)
    shape_r = data["shape_relative"].subs({e: eps*z, k: 1/v})
    shape_s = data["shape_sum"].subs({e: eps*z, k: 1/v})
    f_r = mass*eps*z*sp.sqrt(shape_r)/2
    f_s = mass*sp.sqrt(shape_s)/2
    # H(0)=1 is a normalized linear solution with relative metric O(1/eps).
    assert sp.simplify(sp.limit(eps/f_r, eps, 0)-2*sp.sqrt(v)/(mass*z)) == 0
    # Scaling the linear mode by eps keeps the relative metric amplitude
    # finite, but its *absolute* composite amplitude then tends to zero.
    projection = data["physical_source_relative"].subs({e: eps*z, k: 1/v})/f_s
    assert sp.limit(eps*projection, eps, 0) == 0
    assert sp.simplify(sp.limit(projection, eps, 0)+2*sp.sqrt(v)/mass) == 0


def test_relative_tuning_does_not_cancel_derivatives_or_normalized_tidal_signal():
    data = limiting.general_jets()
    v, a = data["v"], data["a"]
    at = {v: 2, a: 11}
    assert data["theta_relative"].subs(at) == 0
    assert data["N_relative"].subs(at) == sp.Rational(59, 2)
    assert data["mass_squared_prime"].subs(at) == 24
    assert data["mass_squared_second"].subs(at) == -34
    vp, vpp = data["vprime"].subs(at), data["vsecond"].subs(at)
    hsecond = data["heavy_growth_coefficient"].subs(at)
    first = -vp/(2*sp.sqrt(2))
    second = -vpp/(2*sp.sqrt(2))+vp**2/(8*sp.sqrt(2))-sp.sqrt(2)*hsecond
    assert sp.simplify(first+2*sp.sqrt(2)) == 0
    assert sp.simplify(second+10*sp.sqrt(2)) == 0
