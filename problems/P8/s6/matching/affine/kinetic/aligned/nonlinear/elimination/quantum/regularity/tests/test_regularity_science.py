"""Exact higher-derivative, preparation, subtraction and continuum checks."""
import pytest
import sympy as sp
from p8_vector_regularity import estimates, frequency, preparation, spectral, verify


def test_generic_residual_and_projected_readout_algebra():
    verify.affine.certify_residuals({**frequency.algebra_checks(), **spectral.algebra_checks(), **estimates.radial_checks()})


def test_actual_eighth_order_reference_envelopes():
    expected = {"transverse": (5207948696836, 473081396588352),
                "longitudinal": (5165763544817, 474200988057393)}
    for kind, pair in expected.items():
        data = frequency.reference(kind)
        assert (data["residual_over_inverse_frequency_power_upper"],
                data["residual_derivative_over_inverse_frequency_power_upper"]) == pair
        assert all(data["proof_checks"].values())
        assert all(value == 0 for value in data["box_reconstructions"])


def test_unchanged_cutoffs_give_the_stated_three_band_constants():
    data = preparation.constants()["common"]
    assert data["initial_mixing_envelopes"] == {6: 1813229, 8: 906392614, 10: 1536706170025}
    assert data["evolution_mixing_envelope"] == 5347035781757616
    assert all(preparation.proof_checks().values())


def test_band_boundaries_and_large_mass_branch():
    for nu, power in ((1000, 6), (3999, 6), (4000, 8), (10**12-1, 8), (10**12, 10)):
        assert preparation.mixing_bound(1000, nu)["preparation_inverse_frequency_power"] == power
    assert preparation.threshold(10**12) == 8*10**12
    assert preparation.mixing_bound(10**12, 8*10**12)["preparation_inverse_frequency_power"] == 10


def test_first_readout_derivative_independently_matches_physical_weights():
    H, lam, om = spectral.wkb.background()["H"], spectral.wkb.background()["lambda"], spectral.omega2
    for name, weights in spectral.tail.physical_weights().items():
        A, B, d = weights["A"], weights["B"], weights["c1"]-lam/2
        expected = sp.Matrix([[om*(spectral.D0(B)+(2*lam+2*d-3*H)*B),
                               2*om*(B-A), spectral.D0(A)-(2*d+3*H)*A]])
        assert spectral.clean(spectral.rows(name, 1)-expected) == sp.zeros(1, 3)


def test_all_twenty_four_differentiated_subtraction_tails():
    for name in spectral.tail.physical_weights():
        for order in range(6):
            data = spectral.reference_tail(name, order)
            assert data["integrable_tail_certified"]
            assert data["no_lower_Laurent_power"]
            assert data["all_majorant_coefficients_nonnegative"]
    assert all(estimates.proof_checks().values())


def test_lower_reference_failure_is_detected_not_declared_state_divergence():
    assert spectral.low_reference_failure() == sp.Rational(31744, 81)
    assert not spectral.reference_tail("longitudinal_energy", 2, 2)["integrable_tail_certified"]
    assert spectral.reference_tail("longitudinal_energy", 2, 3)["integrable_tail_certified"]


def test_actual_local_derivatives_have_exact_continuous_boxes():
    for observable in ("energy", "pressure"):
        for order in range(6):
            data = estimates.local_derivatives(observable, order)
            assert all(value == 0 for value in data["reconstructions"].values())
    assert estimates.local_derivatives("pressure", 5)["coefficients"][0] == 0


def test_exact_physical_scale_example_and_normalization():
    first = estimates.physical_bounds(10**24, 1000)
    second = estimates.physical_bounds(2*10**24, 1000)
    for name, rows in first["normalized_derivative_bounds"].items():
        for order, row in rows.items():
            assert 0 < row["total"] < sp.Rational(1, 10**18)
            assert row["total"] == 4*second["normalized_derivative_bounds"][name][order]["total"]
    assert first["normalized_derivative_bounds"]["energy"][5]["total"] < sp.Rational(1165, 10**22)
    assert first["normalized_derivative_bounds"]["pressure"][5]["total"] < sp.Rational(1132, 10**22)
    assert not first["perturbed_history_feedback_or_full_quantum_solution_claim"]


@pytest.mark.parametrize("order", (True, 1.0, sp.Float(1), sp.Integer(1)))
def test_native_order_validation_survives_cache_warmup(order):
    spectral.rows("longitudinal_energy", 1)
    with pytest.raises(ValueError, match="native"):
        spectral.rows("longitudinal_energy", order)


def test_all_proof_checks_and_exact_input_controls():
    assert all(verify.proof_checks().values())
    assert verify.controls()["rejected_inputs"] == 120
