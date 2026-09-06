import json
from fractions import Fraction as F

import pytest
import sympy as sp
from p8a_remainder import verify as prior
from p8a_residual import bounds, independent, kernel, plateau


def test_residual_exact_kernel_and_plateau_identities():
    assert all(sp.simplify(value) == 0 for value in kernel.identities().values())
    assert all(sp.simplify(value) == 0 for value in plateau.identities().values())


def test_residual_local_log_and_1922_bounds_have_positive_rational_certificates():
    assert all(value > 0 for value in kernel.elementary_margins().values())
    margins = plateau.elementary_margins()
    assert margins["frozen_ratio_positive_polynomial_coefficients"] == [sp.Rational(1, 32), sp.Rational(97, 992), sp.Rational(2, 961)]
    assert all(value > 0 for value in margins["frozen_ratio_positive_polynomial_coefficients"])
    assert margins["literal_denominator"] == 1922 == 2*31**2
    assert sp.Rational(1, 2048) > sp.Rational(1, 1922**2)


def test_residual_named_length_is_physical_and_kernel_is_actual_clock():
    eta, start, normalization, eta_star = sp.symbols("eta start A eta_star", positive=True)
    a = sp.Function("a", positive=True)(eta)
    length = kernel.named_length(normalization, eta_star)
    assert sp.simplify(length-4*sp.sqrt(2)*(normalization*eta_star**2/2)) == 0
    value = kernel.retarded(a, eta, start, span=3*eta_star, length=length)
    assert value.has(sp.diff(a, eta, 2))
    assert value.has(sp.log(3*a/(2*normalization*eta_star)))


def test_residual_kernel_majorants_have_correct_actual_conformal_dimensions():
    length, delta = sp.symbols("eta_star delta", positive=True)
    b = list(map(sp.Integer, (2, 3, 5, 7)))
    base = kernel.norm_bounds(b, 3, 3, 2, 8)
    actual = kernel.norm_bounds([delta*v/length**(j+2) for j, v in enumerate(b)],
                                3*length, 3, 2/length, 8/length**2)
    assert all(sp.simplify(actual[j]*length**(j+2)/delta-base[j]) == 0 for j in range(3))


@pytest.mark.parametrize("values", [[1, 2, 3], [1, -2, 3, 4], [1.0, 2, 3, 4], [sp.oo]*4])
def test_residual_incomplete_inexact_or_unproved_kernel_norms_rejected(values):
    with pytest.raises(ValueError):
        kernel.norm_bounds(values, 3, 3, 2, 4)


def test_residual_fraction_only_replay_matches_exact_coefficients():
    separate = independent.replay(json.loads(prior.REPORT.read_text()))
    direct = bounds.calibration()
    for name in ("kernel_norm_constants", "effective_S_Born_norm_constants"):
        assert list(map(str, direct[name])) == separate[name]
    for name in ("first_constants", "second_constants"):
        assert {key: str(value) for key, value in direct[name].items()} == separate[name]


def test_residual_rounded_bounds_enclose_all_exact_coefficients():
    data = bounds.calibration()
    for name in bounds.NAMES:
        assert data["first_constants"][name] < data["rounded_first_constants"][name]
        assert data["second_constants"][name] < data["rounded_second_constants"][name]


def test_residual_zero_amplitude_exact_past_and_finite_domain_endpoints():
    assert all(value == 0 for value in bounds.stress_error(0).values())
    assert all(value == 0 for value in bounds.see_residual(0).values())
    cap = bounds.calibration()["delta_bar"]
    assert all(value > 0 for value in bounds.stress_error(cap).values())
    assert all(value > 0 for value in bounds.see_residual(cap).values())


@pytest.mark.parametrize("delta", [-1, 0.1, sp.oo, sp.nan, F(1, 2), sp.Symbol("unknown")])
def test_residual_outside_or_unproved_amplitude_rejected(delta):
    with pytest.raises(ValueError):
        bounds.stress_error(delta)
    with pytest.raises(ValueError):
        bounds.see_residual(delta)


@pytest.mark.parametrize("normalization,eta_star", [(0, 1), (1, 0), (-1, 1), (1, 1.0)])
def test_residual_positive_exact_physical_scales_required(normalization, eta_star):
    with pytest.raises(ValueError):
        bounds.stress_error(0, normalization=normalization, eta_star=eta_star)


def test_residual_exact_units_and_quantum_source_restoration():
    normalization, length, hbar = sp.symbols("A eta_star hbar", positive=True)
    delta = sp.Rational(1, 10**14)
    bare = bounds.stress_error(delta)
    scaled = bounds.stress_error(delta, normalization=normalization, eta_star=length, hbar=hbar)
    for name in bare:
        assert sp.simplify(scaled[name]*normalization**4*length**8/hbar-bare[name]) == 0
    residual_bare = bounds.see_residual(delta)
    residual_scaled = bounds.see_residual(delta, normalization=normalization, eta_star=length)
    for name in residual_bare:
        assert sp.simplify(residual_scaled[name]*normalization**2*length**4-residual_bare[name]) == 0


def test_residual_finite_epsilon_one_example_and_positive_reference_ratios():
    delta = sp.Rational(1, 10**14)
    assert bounds.amplitude(1, delta/16) == delta
    assert all(value < sp.Rational(1, 100) for value in bounds.reference_accuracy(delta).values())
    example = independent.replay(json.loads(prior.REPORT.read_text()))["finite_physical_amplitude_example"]
    assert all(sp.Rational(value) < sp.Rational(1, 10**17)
               for value in example["actual_SEE_residual_over_reference_classical_scales"].values())
    assert not example["positivity_claim_for_arbitrary_Hadamard_states"]
    assert not example["the_metric_is_an_exact_or_nearby_SEE_solution"]


def test_residual_other_numeric_scheme_and_unproved_selector_not_silently_accepted():
    with pytest.raises(TypeError):
        bounds.stress_error(0, gamma=1)
    with pytest.raises(TypeError):
        bounds.stress_error(0, length=1)
    with pytest.raises(TypeError):
        bounds.coefficients(rounded="yes")


def test_residual_exact_bounds_are_no_larger_than_rounded_ones():
    delta = sp.Rational(1, 10**14)
    for function in (bounds.stress_error, bounds.see_residual):
        exact, rounded = function(delta, rounded=False), function(delta, rounded=True)
        assert all(exact[name] < rounded[name] for name in exact)
