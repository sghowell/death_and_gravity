"""Independent Fraction, source-normalization and all-history controls."""

import copy
import json
from fractions import Fraction as F

import pytest
import sympy as sp
from p8a_maxwell_thermal import bounds, dynamics, independent, state, verify


def inputs():
    return (json.loads(verify.prior.REPORT.read_text()),
            json.loads(verify.prior.prior.prior.REPORT.read_text()))


def test_all_calibrations_reconstructed_without_primary_formula_imports():
    result = independent.replay(*inputs())
    for name, primary in (("state", state.calibration()), ("dynamics", dynamics.calibration()),
                          ("bounds", bounds.calibration())):
        assert result[name] == verify.serialize(primary)
    assert not result["imports_new_primary_formula_modules"]


@pytest.mark.parametrize("y,delta", [(F(50, 27), F(1, 10**8)), (F(2), F(1, 10**10)),
                                    (F(3), F(1, 10**12)), (F(100), F(1, 10**8))])
def test_independent_response_polynomial_jets_match_exact_branch(y, delta):
    coupling = F(31, 180)*delta
    expected = independent.jet_values(y, coupling)
    actual = dynamics.branch_point(sp.Rational(y), sp.Rational(delta))
    assert list(map(str, expected)) == list(map(str, (-v for v in actual["normal_Hubble_jets_scaled_by_tau"])))


def test_continuous_response_bounds_rebuilt_from_domain_not_copied_errors():
    total, scheme, lipschitz = independent.response_bounds(F(1, 100))
    assert total == list(map(F, ("16/25", "1728/25", "155136/25", "14770176/25")))
    assert scheme == [64, 6144, 589824]
    assert lipschitz == [8, 96, 1536]
    assert independent.response_bounds(F(1, 200))[0][0] == total[0]/2


@pytest.mark.parametrize("component,key", [("rho", "4,0,0,0"), ("pressure", "2,1,0,0"),
                                           ("EED", "4,0,0,0"), ("trace", "2,1,0,0")])
def test_wrong_sign_or_omitted_physical_anomaly_cannot_enter_reconstruction(component, key):
    _, photon = inputs()
    damaged = copy.deepcopy(photon)
    damaged["derived_constants"]["stress_jet_coefficients"][component]["universal"][key] = "0"
    with pytest.raises(ValueError, match="polynomial"):
        independent.photon_coefficients(damaged)


def test_changed_parent_tube_or_coupling_breaks_actual_history_margin():
    cosmology, photon = inputs()
    damaged = copy.deepcopy(cosmology)
    damaged["derived_constants"]["geometry"]["anchored_C3_slice"]["error_caps"][0] = "0"
    with pytest.raises(ValueError, match="margin"):
        independent.replay(damaged, photon)


def test_independent_high_branch_rejects_positive_algebraic_scale():
    lam = F(31, 18000000000)
    assert F(1, 2) < lam*20000**2 < 1
    with pytest.raises(ValueError, match="low branch"):
        independent.jet_values(20000, lam)
