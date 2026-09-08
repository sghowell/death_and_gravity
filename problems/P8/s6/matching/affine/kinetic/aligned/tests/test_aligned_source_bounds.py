"""Actual coefficient ODE, exact source differentiation and strict domain tests."""
import pytest
import sympy as sp
from p8_affine_aligned import alignment
from p8_affine_aligned import source_bounds as b


def test_original_ODE_and_electric_source_identities():
    assert all(value == 0 for value in b.checks().values())
    assert all(value is True for value in b.proof_checks().values())
    assert b.ode_bridge()["basepoint"] == -1


def test_nonlinear_source_second_variation_independently():
    bg = alignment.parent.old.background()
    n, dk, Q2 = sp.symbols("n delta_K_hat Q_second_order", real=True)
    expected = -2*n*dk/bg["h"]-6*bg["H"]*n**2/bg["h"]+sp.Rational(3, 2)*Q2
    assert sp.factor(alignment.rolling()["second_order"]-expected) == 0
    assert expected.subs({n: 1, dk: 1, Q2: 0, alignment.u: 0}) == -2


def test_exact_source_norm_values_and_zero_controls():
    result = b.bound_values(sp.Rational(1, 10), 2, 3, sp.Rational(1, 100))
    assert result["normal_source_upper"] == sp.Rational(9, 40)
    assert result["electric_source_upper"] == sp.Rational(37, 100)
    assert result["full_action_or_retarded_remainder_claim"] is False
    assert b.bound_values(0, 0, 0, 0)["normal_source_upper"] == 0
    assert b.bound_values(0, 0, 0, 0)["electric_source_upper"] == 0


def test_bound_needs_spatial_clock_control_even_at_a_point_on_the_clock():
    # eta=0 at a point does not set its spatial derivative to zero.
    assert b.bound_values(0, 2, 0, 1)["electric_source_upper"] == 6


@pytest.mark.parametrize("values", ((-1, 0, 0, 0), (sp.Rational(11, 100), 0, 0, 0),
                                   (0, -1, 0, 0), (0, 0, -1, 0), (0, 0, 0, -1),
                                   (0.1, 0, 0, 0), (False, 0, 0, 0), (0, sp.oo, 0, 0)))
def test_source_domain_rejects_bad_or_unproved_norms(values):
    with pytest.raises((TypeError, ValueError)):
        b.bound_values(*values)
