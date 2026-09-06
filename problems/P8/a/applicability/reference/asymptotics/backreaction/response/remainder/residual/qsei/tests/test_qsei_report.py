import copy
import json

import pytest
from p8a_qsei import verify


@pytest.fixture(scope="module")
def replay():
    return verify.build_report()


def test_stored_certificate_is_exactly_replayed(replay):
    verify.validate_report(json.loads(verify.REPORT.read_text()), replay)
    assert replay["claim"] == "P8-A.9"
    assert replay["prior_A8_sha256"] == verify.PRIOR_SHA


def test_exact_source_hash_inventory_includes_independent_audit(replay):
    hashes = replay["source_sha256"]
    assert "tests/test_qsei_independent_audit.py" in hashes
    assert "src/p8a_qsei/independent.py" in hashes
    assert "notes/proof.md" in hashes
    assert "FORMULATION.md" in hashes
    assert all(verify.sha(verify.ROOT/path) == digest for path, digest in hashes.items())


@pytest.mark.parametrize("key", ["auxiliary_flat_kernel_is_a_perturbed_state",
    "alternative_finite_prescription_automatically_covered", "third_sampler_derivative_required",
    "SEE_solution_or_nearby_metric_control", "complete_A1_geodesic_domain_or_initial_curvature_hypothesis",
    "massive_interacting_or_all_realistic_fields", "new_cosmological_incompleteness_or_full_P8a_closure"])
def test_scope_is_not_promoted_beyond_proof(replay, key):
    assert replay["scope"][key] is False


@pytest.mark.parametrize("mutation", ["coefficient", "source_hash", "prior_hash", "state_claim", "identity"])
def test_tampered_report_is_rejected(replay, mutation):
    changed = copy.deepcopy(replay)
    if mutation == "coefficient":
        changed["derived_and_independent_constants"]["all_sampler_constants"]["rounded_QSEI_coefficient"] = "1"
    elif mutation == "source_hash":
        changed["source_sha256"]["notes/proof.md"] = "0"*64
    elif mutation == "prior_hash":
        changed["prior_A8_sha256"] = "0"*64
    elif mutation == "state_claim":
        changed["scope"]["auxiliary_flat_kernel_is_a_perturbed_state"] = True
    else:
        changed["exact_spectral_moments"]["full_UV_moment"] = "1/6"
    with pytest.raises(ValueError):
        verify.validate_report(changed, replay)


def test_changed_prior_certificate_fails_before_replay(monkeypatch):
    original = verify.sha
    monkeypatch.setattr(verify, "sha", lambda path: "0"*64 if path == verify.prior.REPORT else original(path))
    with pytest.raises(ValueError, match="pinned A.8"):
        verify.prior_checks()


def test_build_report_does_not_write_any_file(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("A read-only replay tried to write a file")

    monkeypatch.setattr(verify.Path, "write_text", forbidden)
    monkeypatch.setattr(verify.Path, "write_bytes", forbidden)
    assert verify.build_report()["claim"] == "P8-A.9"
