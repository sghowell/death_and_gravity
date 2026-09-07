import copy
import json

import pytest
from p8_star import verify


@pytest.fixture(scope="module")
def actual():
    return verify.build_report()


def test_star_frozen_report_replays_all_sources_and_ancestors(actual):
    verify.validate_report(json.loads(verify.REPORT.read_text()), actual)


def test_star_immutable_lineage_and_independent_bridges(actual):
    assert actual["prior_context_sha256"]["S6_14_physical_u_and_curvature_lineage"] == verify.PRIOR_SHA
    assert actual["prior_context_sha256"]["S6_13_literal_auxiliary_action"] == verify.ACTION_SHA
    assert actual["prior_context_sha256"]["adopted_S6_conditional_matching_contract"] == verify.CONTRACT_SHA
    bridges = actual["primary_independent_bridges"]
    assert len(bridges["literal_coframe_stresses"]) == 3
    assert [row["leaf_count"] for row in bridges["independent_clock_rate_fixtures"]] == [1, 2, 3, 0]


def test_star_report_retains_algebraic_and_disconnected_exceptions(actual):
    assert "algebraic" in actual["theorem"]["not_global_exact_K"].lower()
    assert "G_u=0" in actual["action_and_domain"]["kinetics"] or "G_u>=0" in actual["action_and_domain"]["kinetics"]
    assert actual["checked_controls"]["actual_disconnected_physical_Hprime_at_zero"] == "2"
    assert actual["checked_controls"]["invalid_domain_calls_rejected"] == 7
    assert "OPEN" in actual["status"]


def test_star_exact_source_inventory_includes_the_reserved_audit_only_once(actual):
    sources = actual["source_sha256"]
    assert "tests/test_star_no_bounce_independent_audit.py" in sources
    assert "notes/proof.md" in sources and "src/p8_star/independent.py" in sources
    assert not any("trimetric/cones/global" in key or key.startswith("../") for key in sources)
    for relative, digest in sources.items():
        assert verify.sha(verify.ROOT/relative) == digest


@pytest.mark.parametrize("omission", ["theorem", "not_established", "prior_context_sha256", "source_sha256"])
def test_star_certificate_whole_claim_or_boundary_omission_is_detected(actual, omission):
    changed = copy.deepcopy(actual)
    del changed[omission]
    with pytest.raises(ValueError):
        verify.validate_report(changed, actual)


def test_star_certificate_one_symbolic_identity_corruption_is_detected(actual):
    changed = copy.deepcopy(actual)
    changed["exact_residuals"]["literal_lapse_scale_variation_and_Bianchi"]["unfactored_bianchi"] = "1"
    with pytest.raises(ValueError):
        verify.validate_report(changed, actual)


def test_star_certificate_source_byte_mutation_is_detected_without_writing(actual):
    changed = copy.deepcopy(actual)
    changed["source_sha256"]["notes/proof.md"] = "0"*64
    with pytest.raises(ValueError):
        verify.validate_report(changed, actual)
