"""Complete ancestry, source manifest and every report-field mutation."""

import copy
import json

import pytest
from p8_offshell_vacuum import verify

KEYS = [
    "schema",
    "claim",
    "date",
    "status",
    "prior_sha256",
    "source_sha256",
    "formulation",
    "written_proofs",
    "exact_residuals",
    "named_exact_check_count",
    "checked_scalar_entries",
    "proof_checks",
    "actual_affine_vacuum_target",
    "off_shell_pointwise_identity",
    "full_centered_resolvent",
    "local_jet_and_substitution_norms",
    "common_Schwartz_domain",
    "exact_calibrations",
    "controls",
    "verdict",
    "not_established",
    "verification_boundary",
]


def test_complete_read_only_replay():
    verify.validate_report(json.loads(verify.REPORT.read_text()), verify.build_report())


def test_own_source_manifest_and_all_report_fields():
    d = verify.build_report()
    assert d["source_sha256"] == {
        str(p.relative_to(verify.ROOT)): verify.sha(p) for p in verify.source_files()
    }
    assert set(d) == set(KEYS)


@pytest.mark.parametrize("key", KEYS)
def test_every_field_mutation_rejected(key):
    d = verify.build_report()
    changed = copy.deepcopy(d)
    changed[key] = "not the certified value"
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, d)


@pytest.mark.parametrize("mutation", ("missing", "extra"))
def test_missing_or_extra_report_field_rejected(mutation):
    d = verify.build_report()
    changed = copy.deepcopy(d)
    if mutation == "missing":
        del changed["common_Schwartz_domain"]
    else:
        changed["full_V_G_B_and_original_P8_closed"] = True
    with pytest.raises(ValueError, match="differs"):
        verify.validate_report(changed, d)


def test_exact_counts_and_original_scope():
    d = verify.build_report()
    assert d["claim"] == "P8-S6.111.OFF_SHELL_VACUUM_QUARTIC_MATCHING"
    assert d["named_exact_check_count"] == d["checked_scalar_entries"] == 24
    assert len(d["proof_checks"]) == 34
    assert d["controls"]["rejected_inputs"] == 35
    assert len(d["source_sha256"]) == 18
