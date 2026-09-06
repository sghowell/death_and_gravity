"""Strict domains, signed sources and conditional-scale controls."""

from fractions import Fraction

import pytest
import sympy as sp
from p8a_maxwell import bounds, focusing
from p8a_maxwell.domain import nonnegative, rational, support_interval


@pytest.mark.parametrize("bad", [True, False, 0.1, sp.Float("0.1"), sp.oo, -sp.oo,
                                 sp.nan, sp.I, sp.Symbol("x"), sp.sqrt(2)])
def test_strict_rational_domain_rejects_nonfinite_or_inexact_input(bad):
    with pytest.raises((TypeError, ValueError)):
        rational(bad)


def test_exact_rationals_and_interval_endpoints():
    assert rational(Fraction(2, 3)) == sp.Rational(2, 3)
    assert support_interval(-2, 3, -1, 1) is True
    for args in ((0, 1, 0, 1), (0, 1, -1, 1), (1, 0, 0, 1)):
        with pytest.raises(ValueError):
            support_interval(*args)
    with pytest.raises(ValueError):
        nonnegative(0, positive=True)
    with pytest.raises(TypeError):
        nonnegative(1, positive=1)


def test_beta_argument_is_mandatory_and_its_sign_remains_visible():
    with pytest.raises(TypeError):
        bounds.envelope(1, 1, 1, 1, 1)
    plus = bounds.envelope(1, 1, 1, 1, 1, beta_m=1)
    minus = bounds.envelope(1, 1, 1, 1, 1, beta_m=-1)
    assert plus["beta_M"] == -minus["beta_M"] == 1
    assert plus["reference_loss_numerator"] == minus["reference_loss_numerator"] == 864
    assert plus["sampler_factor"] == sp.Rational(23, 12)
    assert plus["derivative_coefficient_pi_squared_over_hbar"] == sp.Rational(529, 1152)


def test_physical_source_dictionary_keeps_Lambda_and_extra_negative_EED():
    d = focusing.geometric_constants(0, 0, 0, 0, 1, beta_m=0, kappa=2, hbar=3,
                                     cosmological_constant=5, other_eed_lower=-7)
    assert d["Q2"] == sp.Rational(3, 4)/sp.pi**2
    assert d["Q0"] == 19
    assert d["raw_constant_before_nonnegative_coarsening"] == 19
    improved = focusing.geometric_constants(0, 0, 0, 0, 1, beta_m=0, kappa=2, hbar=3,
                                            cosmological_constant=-5, other_eed_lower=7)
    assert improved["Q0"] == 0
    assert improved["raw_constant_before_nonnegative_coarsening"] == -19


def test_nonnegative_hubble_caps_and_positive_duration_are_required():
    for args in ((-1, 0, 0, 0, 1), (0, 0, 0, 0, 0)):
        with pytest.raises(ValueError):
            bounds.envelope(*args, beta_m=0)


def test_planck_scaling_is_an_identity_not_a_focusing_witness():
    assert all(sp.simplify(value) == 0 for value in bounds.scaling_identities().values())
    d = focusing.calibration()["Hmax_times_tau_at_most_one"]
    assert d["K_trace_times_tau_cap"] < d["A1_zeta_zero_gradient_lower"]
    assert d["small_Q2_implies_focusing"] is False
    assert focusing.calibration()["A1_initial_pointwise_Ricci_premise_removed"] is False


def test_every_classical_source_sign_identity_is_exact():
    assert all(sp.simplify(value) == 0 for value in focusing.identities().values())
