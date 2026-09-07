"""Whole-box derivative bounds, pole cancellation and corner controls."""

from fractions import Fraction as F

import pytest
import sympy as sp
from p8_exact_stationary.intervals import Interval
from p8_own_scattering import potential


def test_entire_box_derivative_and_potential_margins():
    data = potential.calibration()
    assert data["outward_enclosures"]["A_v"] == Interval(F(-128747179, 10**6), F(37954017, 10**6))
    assert data["outward_enclosures"]["pump"] == Interval(F(-21449173, 10**6), F(-2841831, 250000))
    assert all(margin > 0 for margin in data["strict_margins"].values())
    assert data["derived_potential_bound"] == F(345169, 8000)
    assert data["strict_margins"]["potential_below_44"] == F(6831, 8000)
    assert all(potential.checks().values())


def test_polynomial_cancellation_and_positive_remainder_are_exact():
    v, delta = sp.symbols("v delta", nonnegative=True)
    d = 1 + v
    E = (20 + 40*v + 30*v**2 + 8*v**3)/d**4
    Z, D = delta + 8*v, delta + 2*(1 - d**(-4))
    assert sp.cancel(Z - D - v**2*E) == 0
    assert sp.cancel(d**4*(20 - E) - 40*v - 90*v**2 - 72*v**3 - 20*v**4) == 0
    assert sp.cancel(d**4*D - Z - delta*(d**4 - 1) - 12*v**2 - 8*v**3 - 2*v**4) == 0
    controls = potential.pole_polynomials()
    assert controls["Z_minus_D_identity"] == (0, 0, 0, 0, 0, 0)
    assert controls["E_upper"] == (0, 40, 90, 72, 20)


def test_exact_center_A_identity_for_symbolic_c():
    c = sp.Symbol("c", positive=True)
    Q0, b0, z0, L0 = 32/c, sp.Integer(2), 3*(c + 2)**2/(2*c), 6 + 12/c
    bv = -b0*((c - 2)*z0 + L0)/3
    F0 = -2*L0/3
    assert sp.factor(bv + c*(c + 2)) == 0
    assert sp.factor(4*Q0*bv/(F0*b0**2) - 8*c) == 0
    assert sp.factor(2*bv/F0 - c**2/2) == 0


@pytest.mark.parametrize("delta", [F(1, 100), F(1, 625), F(1, 10000)])
def test_exact_positive_delta_center_jets_and_lapse(delta):
    data = potential.center_data(delta)
    c = 2 + delta
    assert data["A"] == 16 + 8*delta
    assert data["b_v"] == -c*(c + 2)
    assert 2*data["b_v"]/data["F"] == c**2/2
    assert all(isinstance(value, F) for value in data.values())


def test_corner_implicit_derivatives_are_exact_not_point_estimates():
    data = potential.center_data(0, extension=True)
    assert (data["zeta_v"], data["zeta_vv"], data["zeta_vvv"]) == (204, 8808, 670824)
    assert (data["b_v"], data["b_vv"], data["b_vvv"]) == (-8, -56, -4800)
    assert (data["F"], data["F_v"], data["F_vv"]) == (-8, -24, -2080)
    assert data["A_v"] == -48
    assert data["ell"] == data["pump"] == -16


def test_direct_kinetic_derivative_reconstructs_the_corner_pump():
    # This contraction uses k=b³/N directly, not potential.py's log formula.
    v = sp.Symbol("v", real=True)
    b = 2 - 8*v - 28*v**2 - 800*v**3
    F0 = -8 - 24*v - 1040*v**2
    N = 2*sp.diff(b, v)/F0
    k = b**3/N
    ell = sp.diff(k, v)/k
    assert sp.cancel(ell.subs(v, 0)) == -16
    assert sp.cancel(sp.diff(ell, v).subs(v, 0)) == -432
    A = 4*(16 - 240*v)*sp.diff(b, v)/(F0*b**2)
    assert sp.cancel(sp.diff(A, v).subs(v, 0)) == -48


def test_u_v_clock_chain_rule_for_the_canonical_pump():
    u, v = sp.symbols("u v", real=True)
    k = 4 - 3*v + 2*v**2
    ell = sp.diff(k, v)/k
    expected = ell + 2*v*sp.diff(ell, v) + v*ell**2
    direct = sp.diff(sp.sqrt(k.subs(v, u**2)), u, 2)/sp.sqrt(k.subs(v, u**2))
    assert sp.cancel(direct - expected.subs(v, u**2)) == 0


def test_two_actual_approach_directions_disprove_joint_continuity():
    v, delta = sp.symbols("v delta", positive=True)
    A_linear = 16 + 8*delta - 48*v
    D_quadratic = delta + 8*v - 20*v**2
    remainder = A_linear/D_quadratic - 16/(delta + 8*v) + 16
    assert sp.limit(remainder.subs(v, 0), delta, 0, dir="+") == 24
    assert sp.limit(remainder.subs(delta, v**2), v, 0, dir="+") == 15
    assert potential.directional_limit(0) == 15
    assert potential.directional_limit(8) == F(73, 4)
    ratio = sp.Symbol("ratio", nonnegative=True)
    formula = 16 + (8*ratio - 48)/(ratio + 8) + 320/(ratio + 8)**2
    assert sp.cancel(formula - (24 - 112/(ratio + 8) + 320/(ratio + 8)**2)) == 0
    assert sp.limit(formula, ratio, sp.oo) == 24


def test_inner_comparison_clock_has_frequency_sqrt7_over2():
    t = sp.Symbol("t", real=True)
    psi = sp.Function("psi")(t)
    y = sp.sqrt(sp.cosh(t))*psi
    residual = (sp.diff(y, t, 2) - sp.tanh(t)*sp.diff(y, t) + 2*y)/sp.sqrt(sp.cosh(t))
    expected = sp.diff(psi, t, 2) + (sp.Rational(7, 4) + sp.Rational(3, 4)/sp.cosh(t)**2)*psi
    assert sp.simplify(residual - expected) == 0


@pytest.mark.parametrize("bad", [True, 0.0, float("nan"), float("inf"), "1/100", None])
def test_exact_domain_and_corner_input_guards(bad):
    with pytest.raises(TypeError):
        potential.center_data(bad)
    with pytest.raises(TypeError):
        potential.directional_limit(bad)


def test_literal_endpoint_and_bad_root_box_rejections():
    with pytest.raises(ValueError):
        potential.center_data(0)
    for delta in (F(-1, 100), F(1, 50)):
        with pytest.raises(ValueError):
            potential.center_data(delta, extension=True)
    with pytest.raises(TypeError):
        potential.center_data(0, extension=1)
    with pytest.raises(ValueError):
        potential.directional_limit(-1)
    with pytest.raises(ValueError):
        potential.partial_jets(root_bracket=Interval(3))
    with pytest.raises(ValueError):
        potential.partial_jets({"v": Interval(0, F(1, 100)), "delta": Interval(0), "zeta": Interval(12)})


def test_scope_does_not_promote_the_bound_to_transfer_or_EFT():
    data = potential.calibration()
    assert not data["combined_remainder_jointly_continuous"]
    assert not data["combined_remainder_derivative_bound"]
    assert not data["literal_delta_zero_allowed"]
    assert not data["fixed_physical_window_transfer_proved"]
    assert not data["low_frequency_EFT_or_UV_verdict"]
