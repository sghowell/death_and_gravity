import sympy as sp
from p8_variable_microlocal import leading, physical, stueckelberg


def test_literal_action_boundaries():
    assert set(stueckelberg.checks().values()) == {0}


def test_stationary_action_algebra():
    assert set(leading.checks().values()) == {0}


def test_vanishing_boundary_coefficient_does_not_allow_omission():
    d = leading.center()
    at = {d["c"]: 4}
    assert d["Fpi_coefficient_at_center"] == 0
    assert d["Fpi_coefficient_derivative_at_center"].subs(at) == 352
    assert d["gradient_without_final_boundary"][0, 0].subs(at) == 448
    assert d["gradient"][0, 0].subs(at) == 800
    assert d["helicity_speed_squared"].subs(at) == sp.Rational(5, 3)


def test_principal_extension_is_not_uniform_time_connection():
    d = leading.center()
    c = d["c"]
    assert sp.factor(d["weight_second_log_jet"]+22+8*c/(c-2)) == 0
    assert sp.limit((c-2)*d["weight_second_log_jet"], c, 2, dir="+") == -16


def test_center_exact_cauchy_map():
    d = physical.equation(0, 4)
    K = stueckelberg.action.K
    expected = 5*(K+60)*(K**2+26*K-24)/(19176*K**3)
    assert sp.factor(d["det_C"]-expected) == 0
    assert d["Mq"].applyfunc(lambda x: sp.limit(x/K, K, sp.oo)) == sp.diag(-sp.Rational(5, 3), -1, -1)
    assert all(physical.degree(x) <= 0 for x in d["Mv"])
    assert d["symplectic"][:3, 3:].applyfunc(lambda x: sp.limit(x, K, sp.oo)) == sp.diag(480, sp.Rational(799, 100), 1)
    assert all(physical.degree(x) <= -2 for x in d["symplectic"][3:, 3:])
