"""Finite massive-log subtraction and continuous complex-domain controls."""
import pytest
import sympy as sp
from p8_vector_bubble import independent, pole, remainder


def test_exact_subtracted_logarithmic_form_and_envelope_proof():
    assert all(value == 0 for value in remainder.checks().values())
    assert all(value is True for value in remainder.proof_checks().values())


def test_scalar_log_remainder_starts_at_six_derivatives():
    z, t = sp.symbols("z t", real=True)
    f0, f2, f4 = sp.symbols("f0 f2 f4", real=True)
    full = (f0+z*f2+z**2*f4)*sp.log(1+t*z)
    after = full-sp.series(full, z, 0, 3).removeO()
    expected = f0*t**3/3-f2*t**2/2+f4*t
    assert sp.simplify(sp.limit(after/z**3, z, 0)-expected) == 0


def test_actual_dimensional_radial_finite_logarithms():
    d = remainder.dimensional_logarithms()
    assert len(d["integrals"]) == 3
    assert len(d["checks"]) == 9
    assert all(value == 0 for value in d["checks"].values())


def test_independent_isotropic_all_order_finite_match():
    assert all(value == 0 for value in independent.isotropic_nonlocal().values())


@pytest.mark.parametrize("order", (3, 4, 5, 8))
def test_independent_isotropic_finite_coefficients_by_direct_parameter_integration(order):
    t = pole.x*(1-pole.x)
    direct = (3*t**order/order-(sp.Rational(1, 2)+6*t)*t**(order-1)/(order-1)
              +(10*t**2-t/2)*t**(order-2)/(order-2))
    reduced = 3*t**order/order-t**(order-1)/(order-1)+t**(order-2)/(4*(order-2))
    assert sp.integrate(direct-reduced, (pole.x, 0, 1)) == 0


def test_feynman_parameter_moment_envelopes_are_exact():
    t = pole.x*(1-pole.x)
    assert sp.integrate(t**2, (pole.x, 0, 1)) == sp.Rational(1, 30)
    assert sp.integrate(t**3, (pole.x, 0, 1)) == sp.Rational(1, 140)
    assert remainder.constants()["integrated_absolute_remainder_constant"] == sp.Rational(4852, 229635)


def test_explicit_mass_scale_and_zero_momentum_controls():
    d = remainder.scale_bound(10**12, 1000, 1)
    assert d["momentum_ball_ratio"] == sp.Rational(1, 10**6)
    assert d["subtracted_loop_kernel_over_reference_density"] == sp.Rational(1, 17280*10**30)
    assert d["full_curved_or_quantum_EFT_remainder"] is False
    assert remainder.scale_bound(10**12, 1000, 0)["subtracted_loop_kernel_over_reference_density"] == 0


@pytest.mark.parametrize("Q", (-1, 2))
def test_outside_complex_momentum_ball_is_rejected(Q):
    with pytest.raises(ValueError):
        remainder.scale_bound(1, 1, Q)
