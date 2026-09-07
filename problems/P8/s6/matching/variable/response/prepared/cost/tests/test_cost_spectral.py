"""Exact arithmetic and omission controls for the analytic cost theorem."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_preparation_cost import spectral


def test_continuous_small_band_and_zero_band():
    assert spectral.lowpass_bounds(100)["tail_coercivity_lower"] == sp.Rational(2, 3)
    assert spectral.lowpass_bounds(0)["exact_trace"] == 0
    assert spectral.lowpass_bounds(0)["tail_coercivity_lower"] == 1
    omega = sp.Symbol("omega", nonnegative=True)
    assert sp.diff(1-omega/300, omega) == -sp.Rational(1, 300)


def test_sinc_remainder_uses_sine_degree_not_an_unproved_alternating_tail():
    result = spectral.sinc_remainder()
    assert result["operator_error_upper"] == sp.Rational(1, 3*39916800)
    assert result["band_gramian_error_upper"] == sp.Rational(16, 25)*result["operator_error_upper"]
    assert result["moment_degrees_required"] == 8
    assert result["actual_band_moments_computed"] is False
    assert spectral.sinc_remainder(0)["operator_error_upper"] == 0
    assert spectral.sinc_remainder(50, 4)["operator_error_upper"] == result["operator_error_upper"]/2**11


def test_fourier_kernel_normalization_and_all_moment_weights():
    x, y, omega = sp.symbols("x y omega", real=True)
    kernel = sp.sin(omega*(x-y))/(sp.pi*(x-y))
    assert sp.limit(kernel, y, x) == omega/sp.pi
    assert len(spectral.moment_weights(100, 4)) == 25
    for left, right, coefficient in spectral.moment_weights(1, 4):
        n = (left+right)//2
        assert sp.simplify(coefficient*sp.pi
                           -(-1)**(n+right)*sp.binomial(2*n, right)/sp.factorial(2*n+1)) == 0


def test_unit_cost_envelope_is_not_a_computed_optimum():
    data = spectral.calibration()
    envelope = data["unit_target_envelope"]
    assert envelope["source_minimum_norm_squared_lower"] == 19**2
    assert envelope["source_minimum_norm_squared_upper"] == sp.Rational(4*10**12, 9)
    assert envelope["spectral_optimizer_norm_squared_upper"] == sp.Rational(2*10**12, 3)
    assert envelope["strict_squared_budget_margin"] == sp.Rational(10**12, 3)
    assert envelope["tail_cost_infimum_lower"] == sp.Rational(722, 3)
    assert data["even_tail_cost_lower"] == 540000
    assert envelope["optimizer_numerically_computed"] is False
    assert envelope["H0_2_attainment_claimed"] is False


def test_minimum_norm_budget_does_not_imply_strict_spectral_budget():
    with pytest.raises(ValueError, match="strict"):
        spectral.minimum_cost_envelope(100, 1, 1, 1, 1)
    with pytest.raises(ValueError, match="strict"):
        spectral.minimum_cost_envelope(0, 1, 1, 1, 1)
    zero = spectral.minimum_cost_envelope(100, 0, 1, 1, 0)
    assert zero["tail_cost_infimum_upper"] == 0


def test_explicit_control_upper_is_not_a_low_frequency_claim():
    bound = spectral.constructive_bounds(100, 1, 2*10**8, 4*10**14)
    assert bound["source_L1_upper"] == 2*10**7
    assert bound["tail_energy_upper"] == 4*10**16
    assert bound["claims_mostly_low_band"] is False
    assert spectral.constructive_bounds(0, 1, 2, 3)["tail_energy_upper"] == 4
    assert spectral.constructive_bounds(100, 1, 2, 3)["tail_energy_upper"] == sp.Rational(9, 10**8)


def test_actual_delta_maps_are_not_replaced_by_their_limit():
    delta = sp.Rational(1, 10**17)
    limit = spectral.finite_delta_error(delta, 1, 10**6)
    prepared = spectral.finite_delta_error(delta, 1, 10**6, target="prepared")
    assert limit-prepared == 200*delta
    assert limit < sp.Rational(1, 1000)
    assert spectral.finite_delta_error(delta, 1, 10**6, sp.Rational(1, 100)) == limit+sp.Rational(21, 50)
    assert spectral.finite_delta_error(sp.Rational(1, 10**9), 1, 10**6) > 1


def test_source_units_use_actual_anisotropic_stress_not_canonical_weights():
    units = spectral.source_units(3, 5, 10)
    assert units["physical_frequency_edge"] == 2
    assert units["physical_L2_squared_multiplier"] == sp.Rational(81, 125)
    assert units["physical_Fourier_energy_multiplier"] == sp.Rational(81, 125)
    assert units["physical_L1_multiplier"] == sp.Rational(9, 5)
    assert units["physical_sup_multiplier"] == sp.Rational(9, 25)


def test_exact_fraction_inputs_and_all_symbolic_checks():
    assert spectral.lowpass_bounds(Fraction(1, 2))["omega"] == sp.Rational(1, 2)
    checks = spectral.checks()
    assert set(checks["residuals"].values()) == {0}
    assert all(value > 0 for value in checks["margins"].values())


@pytest.mark.parametrize("bad", [True, False, 0.1, sp.Float("0.1"), sp.oo, -sp.oo,
                                  sp.nan, sp.Symbol("x", positive=True), "1"])
def test_every_inexact_or_unproved_frequency_is_rejected(bad):
    with pytest.raises(TypeError):
        spectral.lowpass_bounds(bad)


@pytest.mark.parametrize("bad", [-1, 101])
def test_frequency_domain_is_not_silently_extended(bad):
    with pytest.raises(ValueError):
        spectral.lowpass_bounds(bad)


@pytest.mark.parametrize("call", [
    lambda: spectral.sinc_remainder(1, -1),
    lambda: spectral.sinc_remainder(1, sp.Rational(1, 2)),
    lambda: spectral.sinc_remainder(1, 4, -1),
    lambda: spectral.minimum_cost_envelope(1, 1, 0, 1, 3),
    lambda: spectral.minimum_cost_envelope(1, 1, 2, 1, 3),
    lambda: spectral.minimum_cost_envelope(1, 1, 1, 1, 3, 2),
    lambda: spectral.minimum_cost_envelope(1, 0, 1, 1, 3, 1),
    lambda: spectral.constructive_bounds(1, 1, -1, 1),
    lambda: spectral.source_units(0, 1),
    lambda: spectral.source_units(1, 0),
    lambda: spectral.finite_delta_error(0, 1, 1),
    lambda: spectral.finite_delta_error(sp.Rational(1, 10**8), 1, 1),
    lambda: spectral.finite_delta_error(sp.Rational(1, 10**9), 1, 1, momentum_squared=0),
    lambda: spectral.finite_delta_error(sp.Rational(1, 10**9), 1, 1, target="raw_g"),
])
def test_invalid_proof_domains_and_evidence_inputs_are_rejected(call):
    with pytest.raises(ValueError):
        call()
