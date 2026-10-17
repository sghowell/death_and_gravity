"""S1 item 3: P2-0 exact identity chain and D4 (certificates/p2_0_identity.json).

Fast (~1 s).  Every step is an exact sympy identity; each has a negative control that must fail.
"""
import pytest

pytest.importorskip("sympy")

from p2 import p20_identity as P


def test_certificate_replays():
    assert P.verify() == []


@pytest.mark.parametrize("fn", [P.c1_transformed_ode, P.c2_divergence_identity, P.c3_boundary_flux])
def test_identity_and_control(fn):
    assert fn()[0]
    assert not fn(1)[0]


def test_real_imag_elimination():
    ok45, ok_el = P.c4_c5_real_imag_and_elimination()
    assert ok45 and ok_el
    bad45, bad_el = P.c4_c5_real_imag_and_elimination(1)
    assert not bad45 and not bad_el


def test_weight_exponents_zero_mode():
    assert P.c6_weight_positive() and not P.c6_weight_positive(1)
    assert all(P.c7_exponents()) and not any(P.c7_exponents(1))
    assert P.c8_zero_mode_derivative() and not P.c8_zero_mode_derivative(1)


def test_d1prime_machine_checked_parts():
    """C0 (separation of Box_g), C9 (Abel identity), C10 (Abel exponent at the horizons)."""
    assert P.c0_separation() and not P.c0_separation(1)
    assert P.c9_abel() and not P.c9_abel(1)
    assert P.c10_abel_exponent() and not P.c10_abel_exponent(1)
