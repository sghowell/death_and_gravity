"""Read-only certificate for the paired insertion row's finite MS local references."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_fermion_ms_mass import verify as parent

from . import audit, bubble, calibration, ownership, source, tadpole

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-fermion-insertion-ms.json"
PARENT_SHA = "f7dd5bec638832ee8df808daf51ac8efe0cab88f0aafd618ec30b3226dafc0e0"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_fermion_insertion_ms/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen insertion MS parent report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_147_fully_rebuilt": PARENT_SHA,
        "same_GY14_nongravitational_reference_model": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
        "only_older_paired_insertion_local_references_advance": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("An insertion MS proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.148.PAIRED_INSERTION_FINITE_MS_LOCAL_REFERENCES",
        "date": "2026-09-10",
        "status": "OLDER_PAIRED_QUADRATIC_INSERTION_NONLOCAL_AND_FINITE_MS_LOCAL_REFERENCES; NOT_VACUUM_OTHER_FORESTS_FULL_MATCHING_CANONICAL_POLE_HIGHER_LOOPS_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/tadpole.md",
            "notes/bubble.md",
            "notes/forests.md",
            "notes/remainders.md",
            "notes/pole.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "full_finite_MS_local_anchors": serialize(
            {"tadpole": payload(tadpole.data()), "bubble": payload(bubble.data())}
        ),
        "fixed_source_and_forest_ownership": serialize(
            {"source": payload(source.data()), "ownership": payload(ownership.data())}
        ),
        "actual_paired_row_local_enclosure": serialize(payload(calibration.data())),
        "partial_primitive_frontier": serialize(audit.frontier()),
        "controls": serialize(audit.controls()),
        "verdict": "The older scalar_Phi2_W1_F0 row now retains its entire outer MS mass reference after the full inner on-shell fermion forest. The inserted tadpole keeps two exact beta-function anchors and the nonzero-pole finite remainder; the heavy bubble keeps the actual heavy and physical light masses. The fixed assigned H-source counterterm cancels the stationary-heavy local mass term. The resulting combined on-shell mass-reference bound is below 1e-10, while the inherited finite slope and nonlocal bounds are unchanged. Only this row advances. Vacuum rows, other forests, full canonical matching and original P8 remain open.",
        "not_established": [
            "Two vacuum primitive rows and source cross terms in the vacuum reference",
            "Other parameter/counterterm matching insertions and complete canonical two-loop amplitude",
            "Complete two-loop pole/error or higher-loop truncation control",
            "V contour/cut control, finite-gravity G, common-parent B or original P8 closure",
        ],
        "verification_boundary": "Exact beta/Gamma reduction and Laurent finite parts, full regulated two-mass bubble and tadpole remainders, explicit fixed-source cancellation, independent nonzero-regulator integrals and scope mutation controls. Analytic arguments are written proofs, not formalization or independent peer review. Native, ordinary, CLI and direct science retain unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The insertion MS report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.148.PAIRED_INSERTION_FINITE_MS_LOCAL_REFERENCES replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
