"""Fresh ordinary report replay and evidence/scope mutation controls."""
import copy
import json

import pytest
from p8_general_reduction import general, verify


def test_frozen_report_rebuilds_ordinary_ancestry():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_every_source_and_independent_audit_is_pinned():
    actual = verify.build_report()
    sources = {str(path.relative_to(verify.ROOT)): verify.sha(path) for path in verify.source_files()}
    assert actual["source_sha256"] == sources
    assert len(sources) == 14
    assert "tests/test_general_independent_audit.py" in sources
    assert "notes/basis.md" in sources


def test_exact_counts_and_ordinary_domain_controls():
    actual = verify.build_report()
    assert [len(group) for group in actual["exact_residuals"].values()] == [45, 9, 15]
    assert actual["exact_residual_count"] == 69
    assert actual["controls"]["rejected_inputs"] == 22
    assert "ORIGINAL_P8_OPEN" in actual["status"]


@pytest.mark.parametrize("key", ("status", "prior_sha256", "source_sha256", "exact_residuals",
                                  "genuine_formal_S6_deformation", "potential_only_curvature_linear_table",
                                  "boundary_representative_control", "not_established"))
def test_changed_evidence_is_not_accepted(key):
    actual = verify.build_report()
    changed = copy.deepcopy(actual)
    changed[key] = "changed"
    with pytest.raises(ValueError):
        verify.validate_report(changed, actual)


def test_two_distinct_sixth_order_statements_and_degeneracy_control():
    actual = verify.build_report()
    assert not actual["genuine_formal_S6_deformation"]["full_degree_six_action_individually_computed"]
    assert actual["genuine_formal_S6_deformation"]["Einstein_second_variation_unchanged"]
    assert "potential-only" in actual["potential_only_curvature_linear_table"]["scope"]
    assert general.calibration()["double_root_fixture"]["P"] == 0
    assert general.calibration()["negative_P_fixture"]["P"] < 0
    assert actual["controls"]["P_sign_not_a_health_verdict"]
