"""Root-authored coordinate, interval, and limiting-branch scientific audit."""
from fractions import Fraction

import pytest
import sympy as sp
from flint import ctx
from p8_exact_stationary import audit, bridges, intervals
from p8_exact_stationary import stationary as s


def test_independent_four_coordinate_action_and_frozen_profiles():
    assert all(value == 0 for value in bridges.checks().values())
    assert len(bridges.checks()) == 14


def test_literal_einstein_actions_from_coordinate_curvature():
    a = audit.action()
    b, n, q, h = (a[name] for name in ("b", "N", "q", "h"))
    expected_f = -3*audit.M2*b*sp.diff(b, audit.T)**2/n + audit.M2*b**3*sp.diff(h, audit.T)**2/(4*n)
    assert sp.simplify(a["f_ADM"]-expected_f) == 0
    assert sp.simplify(a["g_R"]+sp.diff(q, audit.T)**2/2) == 0
    assert sp.simplify(a["quadratic"].subs({h: q}).doit()
                       -audit.M2*(1+b**3/n)*sp.diff(q, audit.T)**2/4) == 0
    # This last restriction is only an action identity, not controlled locking.


@pytest.mark.parametrize("precision", (128, 256, 384))
def test_independent_arb_encloses_whole_box_and_restores_precision(precision):
    before = ctx.prec
    result = audit.arb_box(precision)
    assert ctx.prec == before
    assert all(value["strict"] for value in result["strict_rational_enclosures"].values())
    assert result["whole_box_not_grid"]
    assert result["endpoint_is_only_analytic_extension"]
    assert result["box"] == intervals.report()["domain"]


@pytest.mark.parametrize("bad", (True, False, 128.0, "256", Fraction(256), 64, 0, -1))
def test_arb_inexact_or_insufficient_precision_is_rejected(bad):
    with pytest.raises((TypeError, ValueError)):
        audit.arb_box(bad)


def test_excluded_action_endpoint_and_regularized_extension_are_distinct():
    p, j = s.profile(), s.joint()
    point = {s.V: 0, s.C: 2, s.ZETA: 12}
    assert p["D"].subs(point) == 0
    assert p["Q"].subs(point) == 16
    assert p["B1"].subs(point) == sp.zoo
    assert sp.simplify(j["N"].subs(point)) == 2
    assert sp.simplify(j["T"].subs(point)) == 12


def test_direct_center_taylor_equations_with_nonzero_lapse_derivative():
    # Substitute independent finite Taylor polynomials into the literal
    # coordinate Euler equations. These exact point audits supplement,
    # rather than replace, the continuous symbolic IFT proof.
    a, f = audit.action(), s.fixed_c()
    for c in (sp.Rational(2001, 1000), sp.Rational(201, 100)):
        time = audit.T
        b = 2+f["b_v_center"].subs(s.C, c)*time**2+f["b_vv_center"].subs(s.C, c)*time**4/2
        lapse = f["N_center"].subs(s.C, c)+f["N_v_center"].subs(s.C, c)*time**2
        profile = s.profile()
        subs = {a["b"]: b, a["N"]: lapse, audit.M2: 1,
                a["p"]: profile["B1"].subs({s.C: c, s.V: time**2}),
                a["w"]: profile["B4"].subs({s.C: c, s.V: time**2})}
        constraint, spatial = (a[name].subs(subs).doit() for name in ("C", "E"))
        for order in (0, 2, 4):
            assert sp.diff(constraint, time, order).subs(time, 0) == 0
        for order in (0, 2):
            assert sp.diff(spatial, time, order).subs(time, 0) == 0


def test_own_f_inner_mass_is_not_the_coupled_relative_mass():
    t = s.tensors()
    e = s.own_equations()
    data = s.inner()
    x = data["x"]
    assert data["own_f_inner_coefficient"] == 16/(1+8*x**2)
    assert data["coupled_relative_inner_coefficient"] == 80/(1+8*x**2)
    assert sp.simplify(t["nu"]*(1/s.M**2+1/t["K_f"])
                       -2*e["beta1"]*(e["b"]+e["N"]/e["b"]**2)/s.M**2) == 0
    assert data["own_f_center_mass_squared"].subs(s.C, 3) == 24/s.TAU**2
    assert data["coupled_algebraic_center_mass_squared"].subs(s.C, 3) == sp.Rational(200, 3)/s.TAU**2
