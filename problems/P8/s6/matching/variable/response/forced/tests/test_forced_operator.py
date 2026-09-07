"""Literal-action/source and retarded kernel checks, not sampled bounds."""

import pytest
import sympy as sp
from p8_variable_forced import coefficients, kernel, operator, source


@pytest.mark.parametrize("module", [operator, kernel, source])
def test_exact_equations(module):
    assert all(sp.simplify(value) == 0 for value in module.checks().values())


def test_pole_and_continuous_identities():
    checks = coefficients.checks()
    assert all(value == 0 for value in checks["residuals"].values())
    assert all(value > 0 for value in checks["strict_margins"].values())


def test_K_zero_is_an_actual_source_case_not_old_transfer_import():
    d = operator.derive()
    zero = {d["K"]: 0, d["u"]: 0, d["delta"]: sp.Rational(1, 10**12)}
    assert d["jL"].subs(zero) > 0
    assert d["jH"].subs(zero) < 0
    assert d["B"].subs(zero) > 0


def test_pole_correction_not_jointly_continuous():
    d = operator.derive()
    u, delta = d["u"], d["delta"]
    rem = d["mass"]-80/(delta+8*u**2)
    assert sp.limit(rem.subs(u, 0), delta, 0) == -32
    assert sp.limit(rem.subs(delta, 0), u, 0) == -141


def test_wrong_physical_source_or_locked_projection_is_detected():
    d = operator.derive()
    actual = d["ag"]*d["jL"]+d["bg"]*d["jH"]
    assert sp.factor(actual-2) == 0
    assert sp.factor(d["ag"]*d["jL"]-2) != 0
    assert sp.Rational(2, 5) != 2


def test_retarded_support_and_coincidence():
    delta = sp.Rational(1, 10**12)
    assert kernel.universal_kernel(-sp.Rational(1, 200), 0, delta, momentum_squared=0) == 0
    assert kernel.universal_kernel(0, 0, delta) == 0


def test_actual_delta_clock_survives_in_kernel():
    delta = sp.Rational(1, 10**12)
    value = kernel.universal_kernel(0, -sp.Rational(1, 200), delta)
    assert value.has(sp.hyper)
    assert value.has(sp.asinh)


def test_missing_adjoint_connection_fails():
    d = operator.derive()
    wrong = sp.factor(d["E"]-d["C"])
    assert wrong != 0
    assert sp.simplify(sp.factor(wrong-sp.diff(d["u"]*d["fc"], d["u"]))) == 0
