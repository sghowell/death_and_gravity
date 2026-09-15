"""Read-only whole original elastic, finite infrared and Proca report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_massive_graviton_cut_isolated_replay import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-massive-elastic-proca-infrared.json"
)
PARENT_SHA = "db850b820bd006cc1400e9e72667ff64ec5d5ea45f1de86c4f1f9683428b9b05"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_massive_elastic_proca_infrared/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen original isolated graviton cut changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_280_complete_isolated_graviton_cut_and_entire_ancestry_rebuilt": PARENT_SHA,
        "S279_remains_rejected_and_archived": True,
        "S275_S276_S277_and_scoped_P8a_unchanged": True,
        "all6_historical_physical_qualifications_unchanged": True,
        "all_original_primitive_and_matching_frontiers_retained": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A massive elastic infrared proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.281.ORIGINAL_COMPLETE_MASSIVE_ELASTIC_ANGULAR_CUT_AND_SPECIFIED_DIMENSIONAL_FINITE_SOFT_MATCHING_WITH_LITERAL_ALL_POLARIZATION_PROCA_SEWING",
        "date": "2026-09-15",
        "status": "SCOPED_FIRST_LOOP_PHYSICAL_ELASTIC_IR_AND_PROCA_CUTS; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/angular.md",
            "notes/dimensional.md",
            "notes/proca.md",
            "notes/sectors.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_elastic_tree": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_finite_soft_Proca_and_species_cut": serialize(
            {
                "packets": {name: payload(packets[name]) for name in names[2:]},
                "all_physical_species": audit.whole_species_cut(),
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The entire original source fixes the complete leading massive elastic tree and full contact/heavy/gravity angular convolution. Its nonforward first-loop finite part is derived with a stated minimal dimensional regulator, full phase-space normalization, evanescent graviton trace and exact matching to the defined S278 analytic soft factor. Literal original Proca stress and all nine physical polarization pairs give the complete nonforward vector cut, calibrated independently to M1. Together with retained M1 and TT kernels these enumerate all physical two-body s cuts below the first heavy threshold without double-counting matter. This is neither the full real amplitude nor crossed forward subtraction, exact finite-G positivity, finite detector or Regge matching, or the original quantum state/domain/bounce. Original V/G/B/P8 remain OPEN.",
        "not_established": [
            "An exact interacting quantum vacuum, LSZ mass or all-orders full-soft-factorization theorem",
            "The complete crossed real amplitude, unphysical low-cut continuation or forward b20",
            "Unknown real local matching, omitted loops, exact detector errors or finite Regge remainder",
            "A unique full off-background dimensional DHOST completion or physical EFT cutoff",
            "Original quantum state/domain/measure/bounce completion or V/G/B/P8 closure",
        ],
        "verification_boundary": "Whole original vacuum jets, exact invariant trees and a complete quadratic angular primitive support the written physical convolution. The dimensional finite part includes an explicit dominated-convergence proof and both phase-space and evanescent contributions, in a named convention. All nine Proca polarization tensors satisfy conservation and completeness, with independent original M1 normalization. Numerical diagnostics are not continuum certificates. Written analysis is not FORMALIZED. Native/direct/ordinary/CLI retain original SymPy; only a captured complete FULL replay uses the separately audited exact-GCD adapter. Frozen science and prior qualifications remain unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete elastic infrared and Proca report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.281 complete original elastic infrared and Proca cut replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
