"""Separate all-time, source and branch audits; no sampled theorem proof."""

import sympy as sp
from p8_bimetric_general import background as bg
from p8_bimetric_monotonic import controls, flow


def test_independent_positive_volume_weight_and_full_euler_bridge():
    t = bg.t
    a, b, ng, nf = bg.a, bg.b, bg.Ng, bg.Nf
    hg, hf = sp.diff(a, t)/(ng*a), sp.diff(b, t)/(nf*b)
    physical_weight = nf*b**3/(ng*a**3)
    grav = -2*bg.MG2*sp.diff(hg, t)/ng-2*bg.MF2*physical_weight*sp.diff(hf, t)/nf
    physical_null = bg.KG/ng**2+physical_weight*bg.KF/nf**2
    source = bg.equations()
    euler = source["EL_a"]-source["EL_Ng"]+physical_weight*(source["EL_b"]-source["EL_Nf"])
    assert sp.factor(euler+grav-physical_null) == 0
    public = flow.derive()
    assert sp.factor(public["full_weighted_null"]-grav) == 0
    assert sp.factor(public["weighted_matter_null"]-physical_null) == 0


def test_off_branch_term_cannot_be_silently_dropped():
    t = sp.Symbol("independent_t", real=True)
    h, y, c, lapse, error = (sp.Function(name)(t) for name in ("h", "y", "c", "N", "e"))
    g, f = sp.symbols("G F", positive=True)
    inertia = g+f*y**2
    hf = h/y+error
    weighted = -2*g*sp.diff(h, t)/lapse-2*f*c*y**3*sp.diff(hf, t)/(c*lapse)
    monotone = -2*inertia**sp.Rational(3, 2)*sp.diff(h/sp.sqrt(inertia), t)/lapse
    assert sp.simplify(weighted-monotone+2*f*y**3*sp.diff(error, t)/lapse) == 0
    assert sp.simplify((weighted-monotone).subs(error, t).doit()+2*f*y**3/lapse) == 0


def test_constant_root_interval_uses_g_equation_only():
    t = sp.Symbol("root_t", real=True)
    h = sp.Function("root_h")(t)
    g, f, y0, lapse = sp.symbols("G F y0 N", positive=True)
    n_g = sp.Symbol("root_n_g", nonnegative=True)
    derivative = sp.diff(h/sp.sqrt(g+f*y0**2), t)/lapse
    actual = derivative.subs(sp.diff(h, t), -lapse*n_g/(2*g))
    assert sp.simplify(actual+n_g/(2*g*sp.sqrt(g+f*y0**2))) == 0
    assert actual.is_nonpositive is True


def test_zero_polynomial_counterexample_is_a_full_decoupled_solution():
    t = sp.Symbol("de_sitter_t", real=True)
    ag, af = sp.exp(t), sp.exp(-t)
    hg, hf = sp.diff(ag, t)/ag, sp.diff(af, t)/af
    assert 3*hg**2-3 == 0
    assert 3*hf**2-3 == 0
    assert 2*sp.diff(hg, t)+3*hg**2-3 == 0
    assert 2*sp.diff(hf, t)+3*hf**2-3 == 0
    z = hg/sp.sqrt(1+(af/ag)**2)
    assert sp.diff(hg, t) == 0
    assert sp.simplify(sp.diff(z, t)).is_positive is True
    assert all(sp.simplify(value) == 0 for value in controls.decoupled_de_sitter_checks().values())


def test_interacting_increasing_H_control_has_actual_canonical_scalar_equation():
    data = controls.increasing_H_example()
    y = data["y"]
    h = sp.sqrt((7*y**2+5/y)/12)
    rho = (7*y**2+5/y)/4+2-15*y/4
    null = -sp.diff(rho, y)/(3*h)
    potential = rho-null/2
    # y'=1 and chi'=sqrt(null)>0 on the stated interval. The field-space
    # derivative of V is V_y/chi', not V_y without the clock conversion.
    scalar_equation = sp.diff(sp.sqrt(null), y)+3*h*sp.sqrt(null)+sp.diff(potential, y)/sp.sqrt(null)
    assert sp.simplify(scalar_equation) == 0
    assert sp.simplify(null-data["n_g"]) == 0
    assert sp.simplify(sp.diff(h, y).subs(y, 1)) == sp.Rational(3, 8)
    assert null.subs(y, 1) == sp.Rational(1, 2)
    assert sp.simplify(sp.diff(h/sp.sqrt(1+y**2), y).subs(y, 1)) == -sp.sqrt(2)/16


def test_degenerate_sign_change_is_detected_away_from_stationary_slice():
    t = sp.Symbol("degenerate_t", real=True)
    g, f = sp.symbols("G F", positive=True)
    h = t**3
    assert sp.diff(h, t).subs(t, 0) == 0
    assert h.subs(t, -1) < 0 < h.subs(t, 1)
    implied_weighted_null = -2*(g+f)**sp.Rational(3, 2)*sp.diff(h/sp.sqrt(g+f), t)
    assert sp.factor(implied_weighted_null+6*(g+f)*t**2) == 0
    assert implied_weighted_null.subs(t, 1).is_negative is True
    # This is a forbidden test profile, not a solution of the parent.


def test_endpoint_only_radius_has_a_sharp_static_control():
    tau = sp.Symbol("tau", positive=True)
    t = sp.Symbol("t", real=True)
    h = 4*t/(tau**2+t**2)
    radius = sp.Rational(8, 5)/tau
    assert sp.simplify(-h.subs(t, -tau/2)-radius) == 0
    assert sp.simplify(h.subs(t, tau/2)-radius) == 0
    assert all(sp.simplify(value) == 0 for value in controls.static_sharpness_checks().values())
    # Error exactly equal to the threshold is attained, not excluded.
