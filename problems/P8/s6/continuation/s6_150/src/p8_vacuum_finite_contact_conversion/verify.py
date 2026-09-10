"""Read-only isolated fixed-contact MS conversion certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_fermion_vacuum import verify as parent

from . import amplitude, audit, calibration, contact, conversion, ownership

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-finite-contact-conversion.json"
PARENT_SHA = "32dcc39d645f6c9375280b99799516dce90c2982b41c694c210bf199793a3909"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_finite_contact_conversion/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen contact-conversion parent report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_149_fully_rebuilt": PARENT_SHA,
        "same_GY14_nongravitational_reference_model": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
        "primitive_rows_unchanged_only_isolated_contact_matching_advances": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A contact-conversion proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.150.ISOLATED_FINITE_CONTACT_MS_CANCELLATION",
        "date": "2026-09-10",
        "status": "ISOLATED_FIXED_CONTACT_NONLOCAL_ORDER_TWO_CANCELLED_WITH_SAME_REGULATED_REFERENCE; NOT_FULL_MS_MATCHING_CANONICAL_POLE_TRUNCATION_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/contact.md",
            "notes/bare.md",
            "notes/amplitude.md",
            "notes/owners.md",
            "notes/bounds.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "exact_regulated_contact_and_bare_map": serialize(
            {
                "contact": payload(contact.data()),
                "conversion": payload(conversion.data()),
            }
        ),
        "full_amplitude_cancellation_and_owners": serialize(
            {
                "amplitude": payload(amplitude.data()),
                "ownership": payload(ownership.data()),
            }
        ),
        "actual_isolated_coordinate_enclosure": serialize(payload(calibration.data())),
        "primitive_and_matching_frontier": serialize(
            {
                "primitive": audit.frontier(),
                "matching": audit.matching(),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The -sigma direction of the fixed finite contact's conversion is discharged at nonlocal order two. The exact affine inverse Lstar+h sigma_D(Lstar)=L and the same entire regulated linear counterterms cancel the full assigned contact insertion against the one-loop re-expansion, retaining both heavy triangles. No extra noncontact G2 or M2 shift or external residue factor is needed for this isolated direction. Its scalar mass insertion is local and fixed by the existing physical mass condition. The exact coordinate inverse is nonsingular at the actual parameters; its tail is not a physical higher-loop bound. Primitive statuses are unchanged. All other matching directions and cross terms, full vacuum/canonical assembly, truncation and V/G/B remain open.",
        "not_established": [
            "Other scale/parameter/field conversion directions or their mixed second-order terms",
            "Complete vacuum/source references, matched two-loop amplitude or canonical pole/error",
            "Physical higher-loop truncation, V contour/cut control or finite-gravity G",
            "Common-parent B or original P8 closure",
        ],
        "verification_boundary": "Exact formal bare-coupling and complete-amplitude expansions with a common holomorphic dimensional contact, independently tested radial integrals, nonconstant affine maps, regulator pole-product mutation and fixed-owner checks. Written analytic arguments and exact replay are not formalization or independent peer review. Native, ordinary, CLI and direct science retain unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The contact-conversion report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.150.ISOLATED_FINITE_CONTACT_MS_CANCELLATION replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
