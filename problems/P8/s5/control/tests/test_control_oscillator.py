import pytest
import sympy as sp
from p8_control import oscillator as o
from p8_s5 import nonlinear_d as m


def test_generic_normal_form_identities():
    assert all(value == 0 for value in o.generic_checks().values())


def test_gamma_bounce_is_finite_and_luminal_at_high_frequency():
    out = o.derive("gamma")
    target = (o.q**3-16*o.q**2+26*o.q+4)/(o.q-2)**2
    assert sp.cancel(out["omega_squared"].subs(o.x, 0)-target) == 0
    assert sp.limit(out["omega_squared"].subs(o.x, 0)/o.q, o.q, sp.oo) == 1
    assert sp.limit(out["mass_correction"].subs(o.x, 0), o.q, sp.oo) == -12


def test_unitary_and_tensor_tail_limits():
    for chart in ("unitary", "tensor"):
        for endpoint in (-1, 1):
            assert o.derive(chart)["omega_squared"].subs(o.x, endpoint) == o.q-6


def test_reciprocal_frequency_derivative_matches_direct_chain_rule():
    f = (1+o.x**3)/(o.q+1)
    for weight in (0, sp.Rational(1, 2), 1, sp.Rational(3, 2)):
        actual = o.derivative(f.subs(o.q, 1/o.r), weight, True)
        expected = o.derivative(f, weight).subs(o.q, 1/o.r)
        assert sp.cancel(actual-expected) == 0


@pytest.mark.parametrize("chart", ["unitary", "gamma"])
def test_against_original_cosmic_time_oscillator_without_compact_derivative(chart):
    old = m.functions()
    u, d = m.u, m.d
    comoving_q = sp.Integer(2500)
    physical_q = comoving_q/d**2
    J, theta, lam = old["J"], old["Theta"], old["Lambda"]
    if chart == "unitary":
        K = J/theta**2
        B = -theta*lam*physical_q/J
        C = 2*(lam**2*physical_q**2/J-physical_q)
    else:
        K = physical_q*J/(physical_q*lam**2-J)
        B = theta*lam*physical_q/J-m.H
        C = 2*theta**2*physical_q**2/J
    f = 3*m.H/2+sp.diff(K, u)/(2*K)
    b = f+B
    point = sp.Rational(3, 4)
    def at(value):
        return value.subs(u, point)
    actual = sp.factor(at(d)*(at(C)/(2*at(K))-at(b)**2-at(sp.diff(b, u))))
    expected = o.derive(chart)["omega_squared"].subs({o.x: sp.Rational(3, 5), o.q: comoving_q/sp.Rational(25, 16)})
    assert sp.cancel(actual-expected) == 0


def test_unproved_denominators_fail_closed():
    with pytest.raises(ValueError, match="Unproved"):
        o.coefficient_bound(1/(1-2*o.x), "tensor")


def test_bound_threshold_meets_every_required_inequality():
    report = o.bounds()
    threshold = sp.Integer(report["q_threshold"])
    assert threshold >= 10**8
    for chart in report["charts"].values():
        C0, C1, C2 = [sp.Rational(v["absolute_bound"]) for v in chart["mass_and_derivative_bounds"]]
        assert threshold >= max(2*C0, C1, C2)
    assert 92/threshold < sp.Rational(1, 100)
    assert 50/threshold < sp.Rational(1, 100)**2


def test_coefficients_obey_majorants_at_exact_points():
    report = o.bounds()
    for chart, point in (("unitary", sp.Rational(1, 8)), ("unitary", 1),
                         ("gamma", 0), ("gamma", sp.Rational(1, 5)), ("tensor", 0)):
        mass = sp.cancel(o.derive(chart)["mass_correction"].subs(o.q, 1/o.r))
        values = (mass, o.derivative(mass, 1, True), o.derivative(o.derivative(mass, 1, True), sp.Rational(3, 2), True))
        for value, entry in zip(values, report["charts"][chart]["mass_and_derivative_bounds"]):
            for inverse_q in (0, sp.Rational(1, 68)):
                evaluated = value.subs({o.x: point, o.r: inverse_q})
                assert abs(evaluated) <= sp.Rational(entry["absolute_bound"])
