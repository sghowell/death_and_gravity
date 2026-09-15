"""Read-only complete heavy-graviton production and soft threshold report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_complete_matter_graviton_endpoint import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT
    / "certificates/polynomial-vacuum-affine-heavy-graviton-production-threshold.json"
)
PARENT_SHA = "191329017c78af9fcb017a3bf77df87bebb2fb47c231948f37dc1e743200c9b7"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(
            ROOT.glob("src/p8_vacuum_affine_heavy_graviton_production_threshold/*.py")
        )
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError(
            "The frozen heavy-graviton production/threshold report changed"
        )
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_290_complete_matter_endpoint_and_entire_ancestry_rebuilt": PARENT_SHA,
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
        raise ValueError("A heavy-graviton production/threshold proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.291.WHOLE_ORIGINAL_HEAVY_GRAVITON_PRODUCTION_FULL_DIMENSIONAL_FORWARD_CUT_AND_SOFT_THRESHOLD_LAURENT_REMAINDER",
        "date": "2026-09-15",
        "status": "SCOPED_HEAVY_GRAVITON_PRODUCTION_AND_THRESHOLD; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/production.md",
            "notes/phase.md",
            "notes/threshold.md",
            "notes/matching.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_complete_production": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_dimensional_cut_and_threshold_Laurent_remainder": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The original four-graph heavy-graviton production tree is Ward complete. Its whole physical TT sew and independent radial/B0 phase normalization determine the full D-dependent forward cut. The positive soft-threshold coefficient, complete first evanescent derivative and compact-window Laurent identity are explicit, with a conditional remainder bound. The unpaired pole is retained; no finite subtraction or virtual cancellation is invented. This is not an exact stable-heavy channel, all mixed four-point cuts, a physical IR/Regge bound or original V/G/B/P8 closure.",
        "not_established": [
            "Paired virtual/pole cancellation or a finite physical threshold observable",
            "Exact stable heavy particle, resonance width or resummation",
            "All mixed four-point cuts/amplitude, full-source b20 or all-loop errors",
            "Finite curved/EFT matching, detector/dressing/IR or Regge remainder",
            "Original state/domain/measure/bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Unchanged source jets and arity, all four tree tensors, full D4/5/6 component projectors plus general-D contraction algebra, independent phase/B0 derivations, complete angular integrals and exact Laurent identities support the written scoped proofs. Numerical angular/evanescent checks are independent cross-checks, not formalized all-domain theorems. Native/direct/ordinary/CLI retain original SymPy; exact-GCD is FULL-only. No frozen source/raw, finite matching, virtual cancellation or width model is changed.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The heavy-graviton production/threshold report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.291 complete heavy-graviton production threshold replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
