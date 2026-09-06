import copy
import json

import pytest
import sympy as sp
from p8a_extended_qsei import focusing, independent, reference, sampling, verify
from p8a_extension import verify as prior


def test_fresh_calibration_and_full_fraction_replay():
    constants = verify.checked_constants()
    duplicate = constants.pop("independent_Fraction_replay")
    assert constants == duplicate
    assert constants["fresh_proper_H2_sampling"]["rounded_absolute_coefficient"] == "2"
    assert constants["actual_comoving_index"]["Q2_over_available_duration_squared_lower"] == "2/5"


def test_all_strict_margins_and_new_span_are_positive():
    for data in (sampling.scattering_calibration(), sampling.calibration(), reference.calibration(), focusing.calibration()):
        assert all(value > 0 for value in data["strict_margins"].values())
    data = sampling.calibration()
    assert data["conformal_envelope"] == sp.Rational(1, 10**6)
    assert data["source_free_conformal_width"] > data["conformal_envelope"]/2
    assert data["difference_coefficient"] < data["derived_absolute_coefficient"] < 2
    assert data["total_spectral_root"] < 1+sp.Rational(2, 10**6)


@pytest.mark.parametrize("value", [-1, sp.Rational(-1, 2)])
def test_absolute_magnitude_rejects_negative_norm_or_hbar(value):
    with pytest.raises(ValueError):
        sampling.absolute_bound(value)
    with pytest.raises(ValueError):
        sampling.absolute_bound(1, hbar=value)


@pytest.mark.parametrize("value", [0.1, sp.Float("0.1"), sp.oo, sp.nan, True])
def test_bound_api_rejects_inexact_nonfinite_or_boolean_inputs(value):
    with pytest.raises((TypeError, ValueError)):
        sampling.absolute_bound(value)


def test_exact_api_amplitude_and_time_dimensions():
    assert sampling.absolute_bound(0) == 0
    assert sampling.absolute_bound(sp.Rational(4, 3), hbar=3) == 1/(2*sp.pi**2)
    assert reference.lower_bound(time_scale=2, hbar=1) == -1/(640*sp.pi**2)
    with pytest.raises(ValueError):
        reference.lower_bound(time_scale=0)


@pytest.mark.parametrize("duration", [0, -1, sp.Rational(4, 10**6)])
def test_index_duration_envelope_guards(duration):
    with pytest.raises(ValueError):
        focusing.index_lower(duration)


def test_index_at_maximum_envelope_is_not_rounded_up():
    span = sp.Rational(3, 10**6)
    lower = focusing.index_lower(span)
    assert lower == 10**6-sp.Rational(3, 8*10**6)
    assert sp.Rational(3, 4) < lower < 10**6


@pytest.mark.parametrize("group", [reference.identities, sampling.identities, focusing.identities])
def test_new_exact_identities(group):
    assert all(sp.simplify(value) == 0 for value in group().values())


@pytest.mark.parametrize("mutation", ["delta", "future_length", "source_rescaled", "distance", "history_derivative"])
def test_independent_replay_rejects_changed_physical_inputs(mutation):
    report = copy.deepcopy(json.loads(prior.REPORT.read_text()))
    data = report["derived_constants"]["complete_weighted_map"]
    if mutation == "delta":
        data["delta"] = "1/100000000000000000000"
    elif mutation == "future_length":
        data["length"] = "1/100000"
    elif mutation == "source_rescaled":
        data["source_off_from"] = "1/2000000"
    elif mutation == "distance":
        data["gate"]["fixed_point_weighted"] = "1"
    else:
        data["geometry"]["history_uprime_cap"] = "1"
    with pytest.raises(ValueError):
        independent.replay(report)
