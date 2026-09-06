import json

import pytest
import sympy as sp
from p8a_qsei import independent, mode_bounds, sampling, verify


@pytest.mark.parametrize("identities", [mode_bounds.identities, sampling.spectral_identities, sampling.clock_identities])
def test_exact_all_sampler_identity_families(identities):
    assert all(value == 0 for value in identities().values())


def test_fraction_polynomial_replay_matches_every_shared_constant():
    separate = independent.replay(json.loads(verify.prior.prior.REPORT.read_text()))
    actual = {**mode_bounds.calibration(), **sampling.calibration()}
    for key, value in actual.items():
        if key in separate:
            assert verify.serialized(value) == separate[key]
    assert separate["IR_Parseval_is_full_complex_not_half"] is True


def test_actual_two_frequency_bound_has_a_uniform_rational_margin():
    data = sampling.calibration()
    assert data["auxiliary_flat_kernel_root_coefficient"] == sp.Rational(2026, 961)
    assert 0 < data["mode_error_root_coefficient_per_delta"] < 64466343
    assert data["maximum_QSEI_coefficient"] < sp.Rational(9, 2) < 5
    assert data["strict_rounding_margin"] > sp.Rational(1, 2)


@pytest.mark.parametrize("delta", [0, sp.Rational(1, 10**16), sampling.DELTA_MAX])
def test_full_closed_amplitude_interval_and_dimensionful_bound(delta):
    assert 0 < sampling.coefficient(delta) < 5
    norm, hbar = sp.symbols("norm hbar", nonnegative=True)
    assert sp.simplify(sampling.absolute_bound(delta, norm, hbar=hbar)
                       -hbar*sampling.coefficient(delta)*norm/(16*sp.pi**2)) == 0
    assert sampling.absolute_bound(delta, 0) == 0


@pytest.mark.parametrize("bad", [-1, 1, 0.0, float("nan"), sp.oo, sp.Symbol("delta", positive=True)])
def test_unproved_or_inexact_amplitude_is_rejected(bad):
    with pytest.raises(ValueError):
        sampling.coefficient(bad)


@pytest.mark.parametrize("bad", [-1, 1.0, sp.oo, sp.Symbol("norm", real=True)])
def test_sampler_norm_requires_a_proved_exact_nonnegative_quantity(bad):
    with pytest.raises(ValueError):
        sampling.absolute_bound(0, bad)


def test_domain_and_scheme_not_silently_expanded():
    with pytest.raises(TypeError):
        sampling.absolute_bound(0, 1, gamma=1)
    with pytest.raises(ValueError):
        independent.replay({"claim": "P8-A.8"})
    with pytest.raises(ValueError):
        verify.zero_checks({"deliberately_wrong_UV_moment": sp.Rational(1, 6)})
