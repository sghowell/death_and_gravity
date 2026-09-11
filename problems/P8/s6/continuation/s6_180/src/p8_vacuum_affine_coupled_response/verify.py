"""Read-only current-parent retained homogeneous tree-plus-vector inverse."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_proca_prepared_inverse import verify as old_curved
from p8_proca_rank_one_inverse import verify as old_scalar
from p8_vacuum_affine_metric_response import verify as parent

from . import audit, classical, inverse, matching, quantum

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-coupled-response.json"
PARENT_SHA = "17a635f2d77b93b8a815cc3e89eae976c06e80d1880c1df02b475d41aec2ad78"
CURVED_SHA = "25eabbd3c7847a746ad7d9b015f27cf63cb3e23aa6d6acfb9b8697fc720bf8fc"
SCALAR_SHA = "0f04fd277040ad21525eb1d9102e6af45c4d695fac3f69d0d37710a7ab8a0db6"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_coupled_response/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError(
            "The frozen current-parent conditional metric-response chain changed"
        )
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    if sha(old_curved.REPORT) != CURVED_SHA or sha(old_scalar.REPORT) != SCALAR_SHA:
        raise ValueError(
            "A source-pinned massive kernel or curved matching proof changed"
        )
    old_curved.validate_report(
        json.loads(old_curved.REPORT.read_text()), old_curved.build_report()
    )
    return {
        "S6_179_fully_rebuilt": PARENT_SHA,
        "S6_87_curved_matching_source_rebuilt": CURVED_SHA,
        "S6_86_massive_scalar_kernel_transitively_rebuilt": SCALAR_SHA,
        "actual_current_tree_and_uncancelled_contacts_rederived_not_old_profile_total_transferred": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A current-parent retained causal inverse gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.180.CURRENT_PARENT_RETAINED_HOMOGENEOUS_CAUSAL_INVERSE_WITH_UNCANCELLED_QUANTUM_CURRENT_AND_BOUNDED_LAPSE_RECOVERY",
        "date": "2026-09-11",
        "status": "UNIQUE_CURRENT_PARENT_PREPARED_RETAINED_LINEAR_CAUSAL_INVERSE_WITH_EXPLICIT_LOCAL_LAPSE_RECOVERY; NOT_SMALL_FULL_INVERSE_QUANTUM_BACKGROUND_STABILITY_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/tree.md",
            "notes/quantum.md",
            "notes/matching.md",
            "notes/inverse.md",
            "notes/recovery.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "actual_current_parent_tree_and_uncancelled_time_channel": serialize(
            {
                "clock": payload(classical.clock_coefficients()),
                "quadratic": payload(classical.quadratic()),
                "tree": payload(classical.adapted()),
                "quantum": payload(quantum.data()),
                "mode_covariance": payload(quantum.mode_covariance()),
            }
        ),
        "matched_retained_inverse_and_physical_reconstruction": serialize(
            {
                "matching": payload(matching.data()),
                "inverse": payload(inverse.data()),
                "first_row": payload(inverse.first_row()),
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching_frontier()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete current affine/CD/M1 target gives a literal prepared matter-reduced tree with a uniformly nonzero fourth-order time channel. Exact time covariance and the full uncancelled conditional vector current leave only one nonlocal scale channel, with no old stress-canceling profile. Whole-dimensional vertices and the shared scheme/state match the existing massive scalar range inverse and curved weak-log remainder. A constructive finite weighted Volterra inverse yields unique smooth prepared retained homogeneous linear responses, and the actual first row gives physical lapse recovery with explicit C1 below15000. The full inverse's kernel majorants remain unevaluated and smallness or stability is not established. The nonzero background residual is not a prepared force, and no self-consistent/full quantum background, cutoff or V/G/B closure follows.",
        "not_established": [
            "A numerical full inverse norm, small backreaction, stability or an all-frequency healthy physical EFT propagator",
            "A compatible quantum preparation, nonlinear corrected background, global-time errors or arbitrary spatial/noise response",
            "Full affine/light/tensor/auxiliary/mixed loops, an interacting measure/state, physical cutoff or omitted-order control",
            "Complete vacuum matching/contour/cuts/truncation, finite-gravity IR/Regge remainder or original V/G/B closure",
        ],
        "verification_boundary": "New literal full-target tree/boundary/charge identities, exact uncancelled time covariance and Ward algebra, full-dimensional current matching, variable-coefficient primitive formulas and continuous local first-row bounds support written matched-kernel and Volterra arguments. Every reused report is rebuilt, with the current tree and contacts independently rederived. Native, direct science, ordinary and CLI use original SymPy; only full regression uses its audited exact GCD adapter. The finite kernel majorants are not numerically evaluated. Not FORMALIZED or a full quantum-parent certificate.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The current-parent retained causal-inverse report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.180.CURRENT_PARENT_RETAINED_HOMOGENEOUS_CAUSAL_INVERSE_WITH_UNCANCELLED_QUANTUM_CURRENT_AND_BOUNDED_LAPSE_RECOVERY replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
