"""Read-only minimal heavy transition and finite formal threshold pair report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_heavy_graviton_production_threshold import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT
    / "certificates/polynomial-vacuum-affine-heavy-resonance-soft-pole-pairing.json"
)
PARENT_SHA = "d95779b399f54eb5b16f11523c9c33713e9c2924f5a221acc2ba88fe77d5b78d"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_heavy_resonance_soft_pole_pairing/*.py"))
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
            "The frozen heavy-graviton production-threshold parent report changed"
        )
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_291_heavy_graviton_production_threshold_and_entire_ancestry_rebuilt": PARENT_SHA,
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
        raise ValueError("A minimal-transition/formal-pairing proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.292.COMPLETE_MINIMAL_HEAVY_TRANSITION_NORMAL_SHEET_MASTERS_AND_FINITE_FORMAL_REAL_VIRTUAL_SOFT_THRESHOLD_PAIRING",
        "date": "2026-09-15",
        "status": "SCOPED_MINIMAL_HEAVY_TRANSITION_AND_FINITE_FORMAL_THRESHOLD_PAIR; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/proper.md",
            "notes/masters.md",
            "notes/pairing.md",
            "notes/matching.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_complete_minimal_transition": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_normal_sheet_masters_and_finite_formal_pairing": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete minimal proper three-point and external-factor graphs have an exact D-dependent master reduction with separated UV and IR origins. Physical-root subtraction fixes the full first evanescent coefficients on the massive normal sheet. The whole Hh real cut and known minimal heavy delta coefficient have a finite formal compact-window pairing with explicit conditional remainder and no chosen finite subtraction. The Coulomb principal value, H-metric/local matching, instability/width, physical IR/Regge and original V/G/B/P8 remain open.",
        "not_established": [
            "Complete H-metric, higher-EFT or local finite matching",
            "Full mixed amplitude, positive physical spectral measure or full-source b20",
            "Exact stable heavy state, width or near-resonance uniformity",
            "Complete physical Coulomb/detector/dressed IR, Regge or all-loop remainder",
            "Original state/domain/measure/bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Original-source graph ownership, full arbitrary-D and component tensor reductions, complete contact/external factors, exact normal-sheet root substitutions, convergent logarithmic derivatives and the full compact-window pairing support the scoped written proofs. Whole-D numerical checks are independent tests, not formalized all-domain theorems. Native/direct/ordinary/CLI use original SymPy; guarded exact-GCD is FULL-only. No frozen source/raw, matching constant, heavy width, or physical IR prescription is changed.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The minimal-transition/formal-pairing report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.292 minimal heavy transition and finite formal threshold pair replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
