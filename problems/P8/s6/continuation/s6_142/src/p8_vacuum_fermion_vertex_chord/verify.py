"""Read-only MS fermion vertex-chord certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_fermion_opposite_chord import verify as parent

from . import audit, calibration, catalog, forest, regions, tail

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-fermion-vertex-chord.json"
PARENT_SHA = "a4b9e871f8c80c7f544d404d9094c1fa668489981c76c5cec57a3e448dff2771"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_fermion_vertex_chord/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen MS fermion vertex-chord report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_141_fully_rebuilt": PARENT_SHA,
        "same_GY14_nongravitational_reference_model": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
        "two_completed_quartic_rows_not_the_whole_sector_or_original_P8": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A MS fermion vertex-chord gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.142.COMPLETE_PAIRED_FERMION_QUARTIC_PRIMITIVES",
        "date": "2026-09-10",
        "status": "TWO_COMPLETE_PAIRED_QUARTIC_PRIMITIVES; NOT_FULL_MATCHED_TWO_LOOP_OTHER_PRIMITIVES_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/catalog.md",
            "notes/forest.md",
            "notes/tail.md",
            "notes/regions.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "proper_MS_vertex_and_region_domains": serialize(
            {"forest": payload(forest.data()), "regions": payload(regions.data())}
        ),
        "primitive_word_ownership": serialize(payload(catalog.data())),
        "soft_tail_and_radial_bound": serialize(payload(tail.data())),
        "actual_completed_quartic_primitive_frontier": serialize(
            {
                "actual": payload(calibration.data()),
                "primitive_frontier": audit.frontier(),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "Every one of the 24 scalar/gauge vertex-chord words is paired with its proper MS Yukawa counterterm. The zero-momentum MS anchor retains all finite dimensional terms. Distinct proved Cauchy radii control the two LOW pieces and the paired HIGH difference. Exact regional integration and complete-subset Lorentz/S4 symmetry bound b2. Adding the disjoint self-energy and opposite-chord groups completes all 60 words in each of two quartic primitive rows. The combined absolute bound is below 1e-1400, relative to the tree below 1e-800. Five other primitive rows, other matching conversions, complete two-loop error, V/G/B and original P8 remain open.",
        "not_established": [
            "Five other wholly unevaluated primitive fermion-sector rows",
            "Other field, parameter and counterterm matching contributions",
            "S6.138 finite outer-reference conversion to common MS",
            "Complete two-loop amplitude, pole or matched error",
            "Higher-loop truncation control, V contours, G, B or original P8 closure",
        ],
        "verification_boundary": "Exact vertex word/cycle and counterterm ownership checks, distinct regional Cauchy domains, explicit convergent radial moments, direct paired-matrix numerical checks and a mutation-checked primitive frontier. Analytic regulator-limit and symmetry arguments are written proofs, not proof-assistant formalization or independent peer review. Native, ordinary, CLI and direct science use unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The MS fermion vertex-chord differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.142.COMPLETE_PAIRED_FERMION_QUARTIC_PRIMITIVES replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
