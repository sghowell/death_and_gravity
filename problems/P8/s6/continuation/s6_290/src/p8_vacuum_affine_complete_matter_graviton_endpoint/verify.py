"""Read-only complete minimal matter-graviton endpoint report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_covariant_gaussian_four_point_matching import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT
    / "certificates/polynomial-vacuum-affine-complete-matter-graviton-endpoint.json"
)
PARENT_SHA = "0a09a6d2e69c5618c6f05f943fc5a8a5da9072e75e897b0a8f430b02e2c98aa1"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_complete_matter_graviton_endpoint/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen covariant Gaussian matching report changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_289_covariant_Gaussian_matching_and_entire_ancestry_rebuilt": PARENT_SHA,
        "S279_remains_rejected_and_archived": True,
        "S275_S276_S277_and_scoped_P8a_unchanged": True,
        "all6_historical_physical_qualifications_unchanged": True,
        "all_original_primitive_and_matching_frontiers_retained": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A complete matter endpoint proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.290.COMPLETE_MINIMAL_MATTER_GRAVITON_ENDPOINT_WITH_UNMATCHED_CURVED_ANCHORS_FULL_PAIR_CUTS_AND_BOUNDED_CROSSED_ENDPOINT_COEFFICIENT",
        "date": "2026-09-15",
        "status": "SCOPED_COMPLETE_MINIMAL_MATTER_ENDPOINT; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/endpoint.md",
            "notes/cuts.md",
            "notes/forward.md",
            "notes/matching.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_complete_matter_endpoint": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_physical_cuts_and_bounded_crossed_endpoint_coefficient": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete original minimal C/g matter-graviton endpoint retains both massive triangles, complete metric/OS/onepoint terms, the full dimensional bubble and explicit curved anchors. Full light/heavy tensor cuts and the whole crossed endpoint coefficient are derived. The known finite-mass graph contribution has magnitude below10^-1005, but the RH residue remains unmatched and other mixed graphs and local coefficients remain. No full-source sign, exact heavy resonance, physical IR/Regge or original V/G/B/P8 closure follows.",
        "not_established": [
            "Finite RH or additional curved/EFT matching coefficients or their bounds",
            "Full mixed matter/gravity four-point amplitude, full-source b20 or all-loop errors",
            "Exact stable heavy atom, physical near-resonance or width resummation",
            "Finite physical detector/dressing/IR and fixed-transfer Regge remainder",
            "Original state/domain/measure/bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact source jets, whole raw/even triangle tensors, dimension-dependent bubble trace, all metric/OS/onepoint terms, complete pair-tensor sewing, literal parameter-cut quadratic, crossed coefficient identities and finite-mass inequalities support the written perturbative proofs. Numerical parameter-versus-sew tests are checks, not formalized all-domain theorems. Native/direct/ordinary/CLI retain original SymPy; exact-GCD is FULL-only. No frozen source/raw or finite curved coefficient is changed.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete matter endpoint report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.290 complete minimal matter endpoint replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
