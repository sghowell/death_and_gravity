"""Exact all-frequency bounds and endpoint/domain exclusion controls."""

import pytest
import sympy as sp
from p8a_see_qsei import scattering


@pytest.mark.parametrize("value", [True, False, 0.1, sp.Float("0.1"), sp.oo, -sp.oo,
                                 sp.nan, "1/0", sp.Symbol("x"), -1])
def test_exact_nonnegative_domain_rejects_unsafe_inputs(value):
    with pytest.raises((TypeError, ValueError)):
        scattering.nonnegative(value)


def test_exact_positive_and_switch_domains():
    assert scattering.nonnegative("3/7") == sp.Rational(3, 7)
    assert scattering.nonnegative(0) == 0
    with pytest.raises(ValueError):
        scattering.nonnegative(0, positive=True)
    with pytest.raises(TypeError):
        scattering.nonnegative(1, positive=1)
    with pytest.raises(ValueError):
        scattering.majorants(0, 0, 0, 0, 0, 0)
    with pytest.raises(ValueError):
        scattering.majorants(2, 1, 1, 1, 1, 1)


def test_volterra_modulus_majorant_and_exact_cap():
    # Termwise cosh(sqrt(B)T)<=exp(BT^2/2), and exp(1/2)<2.
    for n in range(12):
        assert sp.factorial(2*n) >= 2**n*sp.factorial(n)
    data = scattering.majorants(1, 1, 2, 0, 0, 0)
    assert data["modulus_exponent"] == sp.Rational(1, 2)
    assert data["modulus_cap"] == 2
    assert sum(sp.Rational(1, 2)**n/sp.factorial(n) for n in range(10)) < 2


def test_generic_zero_potential_has_zero_scattering_error():
    data = scattering.majorants(3, 0, 0, "1/2", "1/3", "1/3")
    for name in ["q_and_qprime", "A_and_two_derivatives", "forward_two_derivatives", "local_one_derivative"]:
        assert all(value == 0 for value in data[name])
    assert data["qprime_minus_uprime_inverse_k"] == 0
    assert data["infrared_error"] == 0


def test_calibrated_scattering_constants_independently_expanded():
    data = scattering.calibration()
    assert data["potential_cap"] == sp.Rational(3, 10**5)
    assert data["forward_two_derivatives"] == [sp.Rational(27, 200000),
        sp.Rational(3, 40000), sp.Rational(1300081, 20000000000)]
    assert data["local_one_derivative"] == [sp.Rational(9, 400000),
        sp.Rational(500081, 40000000000)]
    assert data["qprime_minus_uprime_inverse_k"] == sp.Rational(3600243, 10**15)
    assert data["infrared_error"] == sp.Rational(117, 200000)
    assert all(value > 0 for value in data["strict_margins"].values())


def test_exact_split_and_nonzero_original_endpoint_control():
    assert all(value == 0 for value in scattering.identities().values())
    k, x = sp.symbols("k x", positive=True)
    q0 = sp.Symbol("q0", nonzero=True)
    direct = q0*(1-sp.exp(-2*sp.I*k*x))/(2*sp.I*k)
    missing_past_endpoint = -q0*sp.exp(-2*sp.I*k*x)/(2*sp.I*k)
    assert sp.simplify(direct-missing_past_endpoint) == q0/(2*sp.I*k)
    assert direct.subs(x, 0) == 0
    assert missing_past_endpoint.subs(x, 0) != 0


def test_naive_backward_supremum_does_not_integrate():
    cutoff = sp.Symbol("cutoff", positive=True)
    assert sp.limit(sp.log(cutoff), cutoff, sp.oo) == sp.oo
    assert scattering.controls()["first_Parseval_is_full_for_complex_products"]


def test_nonlinear_history_remainder_is_at_least_quadratic():
    amplitude = sp.Symbol("amplitude", nonnegative=True)
    t, b, m = sp.symbols("T B M", positive=True)
    l0 = 2*t*amplitude*b
    l1 = amplitude*b+t*amplitude*m+t**2*amplitude**2*b**2
    remainder = sp.expand(amplitude*m*l0+amplitude*b*l1)
    assert remainder.coeff(amplitude, 0) == 0
    assert remainder.coeff(amplitude, 1) == 0
    assert remainder.coeff(amplitude, 2) == b**2+3*b*m*t
    assert remainder.coeff(amplitude, 3) == b**3*t**2
