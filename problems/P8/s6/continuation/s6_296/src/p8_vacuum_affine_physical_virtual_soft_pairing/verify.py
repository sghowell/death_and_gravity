"""Read-only selected physical real/virtual reference-conversion certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-physical-virtual-soft-pairing.json"
)
PARENT_SHA = "c5aeb643e173157f280201af2276924d5a211a5f1f97a271a0ed8dfc8c6d4d5f"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_physical_virtual_soft_pairing/*.py"))
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
            "The frozen complete scalar-bremsstrahlung parent report changed"
        )
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_295_complete_scalar_bremsstrahlung_and_real_error_entire_ancestry_rebuilt": PARENT_SHA,
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
        raise ValueError(
            "A selected physical real-virtual conversion proof gate failed"
        )
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.296.SELECTED_ONE_NEWTON_PHYSICAL_REAL_VIRTUAL_PAIRING_FINITE_ANALYTIC_REFERENCE_CONVERSION_AND_EXPLICIT_ERROR",
        "date": "2026-09-15",
        "status": "SCOPED_SELECTED_PHYSICAL_INCLUSIVE_REFERENCE_CONVERSION_AND_EXPLICIT_ERROR; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/soft.md",
            "notes/continuity.md",
            "notes/inclusive.md",
            "notes/matching.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_dimensional_soft_conversion": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_dimensional_continuity_and_inclusive_reference": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The selected gravitational-dressing real/virtual contribution to the matter-Born rate has a finite conversion to the unchanged S278 analytic hard reference. The full-D angular derivative, phase normalization and correct physical sheet fixDelta_soft=[K1+(EulerGamma-2-lnpi)K0]/(8pi^2 kappa). A fixed-domain positive dimensional continuation supplies the physical nonsoft real remainder, and the total reference-conversion error is below10^-792 at original parameters on the compact domain. Independent hard matching, full gravity, higher loops/operators, amplitude Coulomb/Regge and original V/G/B/P8 remain open.",
        "not_established": [
            "Independent finite hard matching, its positivity or UV completion",
            "Pure-gravity radiation and other hard-loop or higher-EFT sectors",
            "Nonperturbative or full dressed physical S-matrix unitarity",
            "Exchanged regulator/resolution/forward limits or all-energy Regge bound",
            "Original state/domain/measure/bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "The exact full-D current and angular/phase identities, unchanged virtual pole records, positive fixed-domain measure and explicit analytic majorants prove this selected perturbative rate dictionary. Both individual real Fock terms remain divergent; only the paired rate and exact-minus-soft real difference are finite. Independent full-D quadratures calibrate the limits, not rigorous numerical errors. Native/direct/ordinary/CLI use original SymPy; guarded exact-GCD is FULL-only. No finite matching value, frozen source/raw or original closure status is changed.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The selected physical real-virtual reference report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.296 selected physical real-virtual reference conversion replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
