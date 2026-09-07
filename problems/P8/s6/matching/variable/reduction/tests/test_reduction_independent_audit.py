"""Separately authored literal tensor/source and compact-variation audits."""
from fractions import Fraction as Q

import pytest
import sympy as sp
from p8_variable_reduction import tensor


def test_full_coordinate_curvature_and_clock_derivation():
    assert len(tensor.checks()) == 29
    assert set(tensor.checks().values()) == {0}


@pytest.mark.parametrize("index", range(5))
def test_each_clock_invariant_has_no_second_tensor_derivative(index):
    value = tensor.geometry()["Li"][index]
    assert sp.diff(value, sp.diff(tensor.q, tensor.t, 2)) == 0
    assert tensor.euler(value).coeff(sp.diff(tensor.q, tensor.t, 4)) == 0


def test_homogeneous_clock_reparam_does_not_create_fourth_tensor_symbol():
    t, q = tensor.t, tensor.q
    f = sp.Function("f")
    theta = sp.Function("theta")(t)
    density = tensor.derive()["DHOST_tensor_density"].subs(tensor.phi, f(theta)).doit()
    assert tensor.euler(density).coeff(sp.diff(q, t, 4)) == 0


def test_compact_trace_zero_variation_is_not_a_boundary():
    t = tensor.t
    pulse = (1-t**2)**2
    assert all(sp.diff(pulse, t, order).subs(t, end) == 0
               for order in (0, 1) for end in (-1, 1))
    # C1 zero extension is H0^2; smooth approximants preserve this positive integral.
    assert sp.integrate(sp.diff(pulse, t, 2)**2/2, (t, -1, 1)) == sp.Rational(64, 5)


def test_moving_curvature_coefficient_is_not_frozen_before_variation():
    d, t, q = tensor.derive(), tensor.t, tensor.q
    kappa = d["kappa"]
    actual = tensor.euler(kappa*sp.diff(q, t, 2)**2/2)
    expected = (kappa*sp.diff(q, t, 4)+2*sp.diff(kappa, t)*sp.diff(q, t, 3)
                +sp.diff(kappa, t, 2)*sp.diff(q, t, 2))
    assert sp.expand(actual-expected) == 0
    assert sp.expand(actual-kappa*sp.diff(q, t, 4)) != 0


@pytest.mark.parametrize("c", (Q(2000000001, 10**9), Q(5, 2), Q(3), Q(4)))
def test_fraction_center_spring_and_curvature_calibration(c):
    b1 = Q(32)/(c*(c-2))
    relative_spring = 2*b1*2
    hidden_Einstein_coefficient = Q(4)
    kappa_from_schur = hidden_Einstein_coefficient**2/(2*relative_spring)
    assert kappa_from_schur == c*(c-2)/16 > 0
    assert kappa_from_schur/2 == c*(c-2)/32


def test_source_preserving_schur_coefficient_includes_physical_g_weight():
    G, F, nu, D = sp.symbols("G F nu D", positive=True)
    kernel = sp.Matrix([[G*D+nu, -nu], [-nu, F*D+nu]])
    actual = sp.factor(kernel[0, 0]-kernel[0, 1]**2/kernel[1, 1])
    local = (G+F)*D-F**2*D**2/nu
    remainder = F**3*D**3/(nu*(nu+F*D))
    assert sp.cancel(actual-local-remainder) == 0
    massive = nu*(G+F)/(G*F)
    response = 1/((G+F)*D)+F/(G*(G+F)*(D+massive))
    assert sp.cancel(kernel.inv()[0, 0]-response) == 0
    # The curvature combination supplies twice kappa in this kernel convention.
    kappa = F**2/(2*nu)
    assert sp.expand(local-(G+F)*D+2*kappa*D**2) == 0


def test_inverse_metric_minimal_matter_variation_has_positive_half_T_sign():
    epsilon = sp.Symbol("epsilon", real=True)
    Y, ww, wv = sp.symbols("Y ww wv", real=True)
    # Rank-one inverse-metric variation Z=w w^T. The determinant lemma is exact.
    density = (Y+epsilon*wv**2)/(2*sp.sqrt(1+epsilon*ww))
    stress_contraction = wv**2-Y*ww/2
    assert sp.diff(density, epsilon).subs(epsilon, 0) == stress_contraction/2
    assert sp.diff(density, epsilon).subs(epsilon, 0) != -stress_contraction/2


def test_curvature_redefinition_induces_nonfree_matter_contact():
    # Contractions of T_chi = v v - g Y/2 in four dimensions.
    Y, Rvv, R, kappa, Planck2 = sp.symbols("Y Rvv R kappa Planck2", real=True)
    T_Ricci = Rvv-R*Y/2
    T_trace = -Y
    induced = kappa/Planck2*(T_Ricci-R*T_trace/6)
    assert sp.expand(induced-kappa/Planck2*(Rvv-R*Y/3)) == 0
    # Pure-chi part of the leading equation Q G=T gives Ricci=v v/Q.
    on_shell_chi = induced.subs({Rvv: Y**2/Planck2, R: Y/Planck2})
    assert sp.factor(on_shell_chi) == 2*kappa*Y**2/(3*Planck2**2)
