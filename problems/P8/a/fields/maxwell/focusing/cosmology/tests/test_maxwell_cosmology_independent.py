"""Independent pinned-polynomial and exact arithmetic omission controls."""

import copy
import json

import pytest
from p8a_maxwell_cosmology import independent, verify


def inputs():
    return (json.loads(verify.prior.REPORT.read_text()),
            json.loads(verify.prior.prior.REPORT.read_text()))


def test_all_new_primary_datasets_match_fraction_only_reconstruction():
    data = verify.checked_constants()
    assert not data["independent_Fraction_replay"]["imports_new_primary_formula_modules"]
    assert data["calibration"]["universal_gate"]["margin_above_one_eighth"] == "3329/350000"


def test_omitting_a_photon_anomaly_monomial_changes_the_independent_cost():
    theorem, photon = inputs()
    expected = independent.replay(theorem, photon)
    damaged = copy.deepcopy(photon)
    damaged["derived_constants"]["stress_jet_coefficients"]["EED"]["universal"].pop("2,1,0,0")
    actual = independent.replay(theorem, damaged)
    assert actual["calibration"]["expanded_affine_cost"]["C0"] != expected["calibration"]["expanded_affine_cost"]["C0"]


def test_finite_beta_reference_term_is_not_silently_dropped():
    theorem, photon = inputs()
    expected = independent.replay(theorem, photon)
    damaged = copy.deepcopy(photon)
    damaged["derived_constants"]["stress_jet_coefficients"]["EED"]["beta_M"].pop("0,0,0,1")
    actual = independent.replay(theorem, damaged)
    assert actual["calibration"]["expanded_affine_cost"]["Cbeta"] != expected["calibration"]["expanded_affine_cost"]["Cbeta"]


def test_tampered_cubic_or_wrong_geometric_weight_is_rejected():
    theorem, photon = inputs()
    bad_moment = copy.deepcopy(theorem)
    bad_moment["derived_constants"]["cubic_moments"]["zeroth"] = "1/3"
    with pytest.raises(ValueError, match="cubic"):
        independent.replay(bad_moment, photon)
    bad_weight = copy.deepcopy(photon)
    bad_weight["derived_constants"]["stress_jet_coefficients"]["EED"]["universal"]["3,0,0,0"] = "1"
    with pytest.raises(ValueError, match="weight"):
        independent.replay(theorem, bad_weight)
