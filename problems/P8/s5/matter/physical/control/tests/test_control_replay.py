import copy
import json

import pytest
from p8_m1_control import audit, verify


@pytest.mark.parametrize("chart,u_value", audit.POINTS)
def test_full_main_oscillator_matches_original_cosmic_time(chart, u_value):
    actual = audit.point_check(chart, u_value)
    assert set(actual["exact_residuals"]) == {
        "full_finite_q_potential", "antisymmetric_connection",
        "canonical_momentum_boundary", "physical_kinetic_whitening"}
    assert all(entry == "0" for rows in actual["exact_residuals"].values()
               for row in rows for entry in row)


def test_symbolic_bounce_is_not_only_a_high_frequency_limit():
    actual = audit.bounce_check()
    assert len(actual["exact_residuals"]) == 19
    assert set(actual["exact_residuals"].values()) == {"0"}


def test_read_only_exact_certificate_replay():
    expected = json.loads(verify.REPORT.read_text())
    verify.validate_report(expected, verify.build_report())
    assert expected["uniform_exact_majorants_and_free_energy_contract"]["q_threshold"] == "100000000000000000000"


def test_certificate_scope_and_energy_mutations_rejected():
    actual = verify.build_report()
    altered = copy.deepcopy(actual)
    altered["not_established"] = []
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(altered, actual)
    altered = copy.deepcopy(actual)
    altered["uniform_exact_majorants_and_free_energy_contract"]["full_window_energy_ratio"] = ["1", "1"]
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(altered, actual)


def test_changed_prior_certificate_rejected(monkeypatch):
    original = verify.sha
    monkeypatch.setattr(verify, "sha", lambda path: "wrong" if path == verify.prior.REPORT else original(path))
    with pytest.raises(ValueError, match="certificate changed"):
        verify.prior_checks()
