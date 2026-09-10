"""Read-only complete paired mixed quartic primitive certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_fermion_outer_ms import verify as parent

from . import audit, bounds, catalog, dimensional, joint, reference

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-fermion-mixed-quartic.json"
PARENT_SHA = "4ec6b42374a9d0bf3832d32389cd73dcb3a9e0de4ce1cdd1fbdd2a52f7215026"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_fermion_mixed_quartic/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen mixed-quartic parent report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_143_fully_rebuilt": PARENT_SHA,
        "same_GY14_nongravitational_reference_model": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
        "only_one_mixed_quartic_primitive_advances": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A mixed-quartic proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.144.COMPLETE_PAIRED_MIXED_FERMION_QUARTIC",
        "date": "2026-09-10",
        "status": "ALL_FOUR_QUARTIC_PRIMITIVE_ROWS_BOUNDED; NOT_FULL_MATCHED_TWO_LOOP_OTHER_CONVERSIONS_HIGHER_LOOP_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/catalog.md",
            "notes/joint.md",
            "notes/zero_soft.md",
            "notes/finite_ms.md",
            "notes/conversion.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_mixed_catalog": serialize(payload(catalog.data())),
        "joint_soft_remainder": serialize(payload(joint.data())),
        "zero_soft_MS_reference_and_conversion": serialize(
            {
                "dimensional": payload(dimensional.data()),
                "finite": payload(reference.data()),
            }
        ),
        "actual_enclosure_and_frontier": serialize(
            {"actual": payload(bounds.data()), "primitive_frontier": audit.frontier()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The full mixed fermion box/tree-Hessian primitive is split exactly into a nonzero-soft remainder and a zero-soft reference. Joint two-loop Schwinger bounds control the remainder without assigning a large Cauchy radius to light lines. The zero-soft kernel is a fixed-mu second mass derivative with a signed spectral density. Exact regulated soft and hard reference integrals retain finite pole products; the nonzero mass-ratio correction is bounded. A literal parent map includes the finite outer MS conversion. All four quartic primitive rows now have bounds, but other forests, matching and canonical contributions, the two-loop pole/error and original P8 remain open.",
        "not_established": [
            "Four remaining primitive vacuum/quadratic rows",
            "Other finite field, parameter and counterterm matching contributions",
            "Complete matched two-loop amplitude, canonical pole or error",
            "Higher-loop truncation control, V contours, finite-gravity G, common-parent B or original P8 closure",
        ],
        "verification_boundary": "Exact operator derivatives, complex denominator gaps, Schwinger/Gamma and Laurent identities, independent numerical integrals and mutation-checked scope. Analytic convergence and continuation arguments are written proofs, not proof-assistant formalization or independent peer review. Native, ordinary, CLI and direct science retain unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The mixed-quartic report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.144.COMPLETE_PAIRED_MIXED_FERMION_QUARTIC replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
