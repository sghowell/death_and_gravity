"""Root-owned independent action, physical-clock, ODE and hierarchy audit."""

from fractions import Fraction as Q

import pytest
import sympy as sp
from p8_composite_light import background, bounds, canonical
from p8_composite_modes import model as m


def test_literal_tensor_action_uses_physical_proper_time():
    data = canonical.physical_coefficients()
    kg = m.G*m.a**3*m.s/8
    kf = m.F*m.a**3*m.y**3*m.s/(8*m.c)
    assert sp.cancel(data["K1"]-kg) == 0
    assert sp.cancel(data["K2"]-kf) == 0
    # dt=dT/Ne and d/dt=Ne*d/dT: multiplying the kinetic, but
    # dividing the stiffness, is essential to the same-clock gap.
    mu = m.y*(m.M4*(m.betas[1]+m.betas[2]*(m.y+m.c)+m.betas[3]*m.c*m.y)
              -m.alpha*m.beta*m.r*m.s*m.p)
    stiffness = m.a**3*mu/(8*m.s)
    relative = kg*kf/(kg+kf)
    assert sp.cancel(data["relative_stiffness_T"]-stiffness) == 0
    assert sp.cancel(data["mass_squared"]-stiffness/relative) == 0


def test_direct_rotating_weight_transformation_and_spatial_cross_term():
    k, w = sp.symbols("K w", positive=True)
    vd, dd, wd, delta = sp.symbols("vd dd wd delta", real=True)
    v, cg, cf = sp.symbols("v cg2 cf2", real=True)
    hd, gd = vd-w*dd-wd*delta, vd+(1-w)*dd-wd*delta
    literal_kinetic = k*(1-w)*hd**2+k*w*gd**2
    assert sp.expand(literal_kinetic-k*(vd-wd*delta)**2-k*w*(1-w)*dd**2) == 0
    literal_gradient = k*(1-w)*cg*(v-w*delta)**2+k*w*cf*(v+(1-w)*delta)**2
    cross_coefficient = sp.expand(literal_gradient).coeff(v, 1).coeff(delta, 1)
    assert sp.factor(cross_coefficient-2*k*w*(1-w)*(cf-cg)) == 0
    # Freezing w before the coordinate change deletes a genuine coupling.
    omitted = sp.expand(literal_kinetic-k*vd**2-k*w*(1-w)*dd**2)
    assert omitted != 0


def test_canonical_normalization_and_boundary_from_independent_functions():
    t = sp.Symbol("t", real=True)
    fs, fr, w = (sp.Function(n, positive=True)(t) for n in ("fs", "fr", "w"))
    ll, hh = sp.Function("ell")(t), sp.Function("heavy")(t)
    # fR/fS=sqrt(w1*w2). Keep wdot instead of a frozen eigenvector.
    literal = fs**2/2*(sp.diff(ll/fs, t)-sp.diff(w, t)*hh/fr)**2
    literal += fr**2/2*sp.diff(hh/fr, t)**2
    om = fs*sp.diff(w, t)/(2*fr)
    normalized = (sp.diff(ll, t)-sp.diff(fs, t)*ll/fs-2*om*hh)**2/2
    normalized += (sp.diff(hh, t)-sp.diff(fr, t)*hh/fr)**2/2
    assert sp.simplify(literal-normalized) == 0
    data = canonical.generic_action()
    boundary = sp.diff(data["normalization_boundary"], data["t"])
    assert sp.expand(data["pre_boundary_action"]-data["post_boundary_action"]-boundary) == 0


def test_adjoint_surface_term_and_retarded_heavy_data_cannot_be_dropped():
    t = sp.Symbol("t", real=True)
    f, g, om, ts, cross = (sp.Function(n)(t) for n in ("f", "g", "om", "ts", "cross"))
    left = f*canonical.B(g, t, om, ts, cross)
    right = canonical.B_adjoint(f, t, om, ts, cross)*g
    assert sp.expand(left-right-sp.diff(2*om*f*g, t)) == 0
    homogeneous = sp.cos(3*t)
    assert sp.diff(homogeneous, t, 2)+9*homogeneous == 0
    omitted_force = canonical.B_adjoint(homogeneous, t, sp.S.One, sp.S.Zero, sp.S(2))
    assert sp.expand(omitted_force-6*sp.sin(3*t)-2*sp.cos(3*t)) == 0
    assert omitted_force.subs(t, 0) == 2


def test_inverse_mass_representative_retains_full_differential_residual():
    t = sp.Symbol("t", real=True)
    light = t**3+t
    mass, extra = 2+t**2, 1-t
    force = 2*(1+t)*(sp.diff(light, t)-t*light)+t**2*light
    representative = canonical.mass_led_pair(light, t, 1+t, t, t**2, mass)
    assert sp.cancel(representative+force/mass) == 0
    residual = sp.diff(representative, t, 2)+(mass+extra)*representative+force
    assert sp.cancel(residual+sp.diff(force/mass, t, 2)+extra*force/mass) == 0
    assert residual.subs(t, 0) != 0
    with pytest.raises(ValueError, match="nonzero"):
        canonical.mass_led_pair(light, t, 1+t, t, t**2, 0)
    assert canonical.generic_action()["mass_squared"].is_positive is not True


def test_full_undivided_reconstruction_generates_symmetric_second_jets():
    y, rho, h = sp.symbols("yy density hubble", positive=True)
    acceleration = background.A
    r = 1+y
    gx = sp.sqrt(y**2+r**3*rho/3)
    fy = -sp.sqrt(y**-2+r**3*rho/(3*y**3))
    ratio = (gx-r*h)/(r*h-y*fy)
    yd = y*r*(gx*fy-h*(gx+fy))/(gx-y*fy)
    rd = -3*h*(rho+2*y/r**2)
    field = (yd, rd, acceleration)
    variables = (y, rho, h)

    def derivative(value):
        return sum(sp.diff(value, x)*v for x, v in zip(variables, field, strict=True))

    at = {y: 1, rho: sp.Rational(1, 2), h: 0}
    root = sp.sqrt(sp.Rational(7, 3))
    results = {"y_prime": yd, "y_second": derivative(yd),
               "c_prime": derivative(ratio), "c_second": derivative(derivative(ratio)),
               "rho_second": derivative(rd)}
    expected = {"y_prime": -root, "y_second": sp.Rational(7, 3),
                "c_prime": -(5+12*acceleration)/(3*root),
                "c_second": (5+12*acceleration)**2/21,
                "rho_second": -3*acceleration}
    public = background.initial_jets(1)
    for name, expression in results.items():
        direct = sp.simplify(expression.subs(at))
        assert sp.simplify(direct-expected[name]) == 0
        assert sp.simplify(public[name]-direct) == 0


def test_zero_algebraic_stiffness_is_not_a_zero_routh_frequency():
    data = background.initial_jets(1)
    at = {background.A: sp.Rational(4, 25)}
    assert data["mass_squared"] == 0
    assert data["Routh_frequency_squared"].subs(at) == sp.Rational(153, 100)
    assert data["heavy_diagonal_at_zero_k"].subs(at) < 0
    assert data["omega"].subs(at) != 0
    # This is a fixed-charge q=0 countercontrol, not a finite-band gap proof.
    assert sp.simplify(data["Routh_frequency_squared"]+data["N_relative"]) == 0


def test_exact_asymmetric_mass_and_all_A_mixing_screen():
    root_lo, root_hi = Q(2102379, 100000), Q(2102380, 100000)
    assert root_lo**2 < 442 < root_hi**2
    mass_lo, mass_hi = (107-5*root_hi)/42, (107-5*root_lo)/42
    assert 0 < mass_lo < mass_hi < Q(9, 200)
    # Uniform on 0<=A<=1/4, using the deliberately looser sqrt442<22.
    numerator = Q(1037)-(376+16*Q(22))/4
    omega2_floor = numerator**2/(16*(60996+3047*Q(22)))
    assert numerator > 0 and omega2_floor > Q(1, 3)
    # On A>=1/4 the physical Hdot/m²=A itself supplies the floor.
    assert Q(1, 4)/Q(9, 200) == Q(50, 9)
    actual = background.initial_jets(2)
    assert sp.simplify(actual["mass_squared"]-(107-5*sp.sqrt(442))/42) == 0
    assert bounds.build()["strict_mass_led_ratio_lower"] == Q(50, 9)


def test_normalization_generated_frequency_has_no_hierarchy_at_these_jets():
    data = background.initial_jets(1)
    A = background.A
    difference = 3*data["N_sum"]-data["Routh_frequency_squared"]
    assert sp.simplify(difference-(36*A**2+30*A+1)/7) == 0
    assert sp.simplify(data["N_sum"]-sp.Rational(12, 7)*(A-sp.Rational(1, 48))**2-sp.Rational(51, 64)) == 0
    # A second independent, tighter rational sqrt442 enclosure checks every
    # coefficient sign in the asymmetric all-A polynomial argument.
    root = sp.Symbol("radical")
    lo, hi = Q(2102379, 100000), Q(2102380, 100000)
    other = background.initial_jets(2)
    for name, sign, constant in (("N_sum", 1, Q(1, 2)), ("Routh_frequency_squared", -1, Q(3, 2))):
        poly = sp.Poly(other[name], A)
        assert poly.degree() == 2
        for j in range(3):
            aff = sp.Poly(sp.radsimp(poly.nth(j)).subs(sp.sqrt(442), root), root)
            assert aff.degree() <= 1
            b, c = (Q(str(aff.nth(i))) for i in range(2))
            ends = [b+c*lo, b+c*hi]
            if j == 0:
                ends = [v-constant for v in ends]
            assert min(sign*v for v in ends) > 0
