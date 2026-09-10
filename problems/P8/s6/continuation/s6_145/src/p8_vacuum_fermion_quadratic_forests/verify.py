"""Read-only paired nonlocal quadratic primitive certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_fermion_mixed_quartic import verify as parent

from . import audit, calibration, catalog, forest, pole, regions

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-fermion-quadratic-forests.json"
PARENT_SHA = "32208fdf95cbaa1bfdef6e3d43f84450b0b2d41e9ec043fc544b6eeca3e5ea0d"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_fermion_quadratic_forests/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen quadratic parent report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_144_fully_rebuilt": PARENT_SHA,
        "same_GY14_nongravitational_reference_model": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
        "two_nonlocal_quadratic_rows_not_finite_local_references": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A quadratic-forest gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.145.PAIRED_QUADRATIC_NONLOCAL_REMAINDERS",
        "date": "2026-09-10",
        "status": "BOTH_QUADRATIC_PRIMITIVE_NONLOCAL_OS_REMAINDERS_BOUNDED; NOT_FINITE_MS_LOCAL_REFERENCES_FULL_MATCHING_CANONICAL_POLE_HIGHER_LOOPS_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/catalog.md",
            "notes/forest.md",
            "notes/regions.md",
            "notes/projection.md",
            "notes/bounds.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_proper_cycle_catalog": serialize(payload(catalog.data())),
        "paired_MS_forest_and_regions": serialize(
            {"forest": payload(forest.data()), "regions": payload(regions.data())}
        ),
        "on_shell_nonlocal_projection": serialize(payload(pole.data())),
        "actual_enclosure_and_partial_frontier": serialize(
            {
                "actual": payload(calibration.data()),
                "primitive_frontier": audit.frontier(),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "Each scalar/gauge quadratic primitive has two self-energy words and one vertex word. The complete proper-cycle catalog includes the whole four-fermion cycle; its local term and the overall affine polynomial vanish under the nonlocal on-shell projection. Both overlapping Yukawa vertex subtractions and both finite MS anchors are retained. Symmetric middle/high momentum regions use distinct legitimate radii for the paired kernel and its complementary subtraction. Exact radial moments and a second Cauchy estimate bound the divided on-shell nonlocal remainder below 1e-799. Finite MS local references, other matching/canonical contributions, full two-loop pole/error and original P8 remain open.",
        "not_established": [
            "Finite MS mass and slope references for these quadratic rows",
            "Two primitive vacuum rows and other counterterm/matching insertions",
            "Complete matched/canonical two-loop amplitude, pole or error",
            "Higher-loop truncation control, V contours, finite-gravity G, common-parent B or original P8 closure",
        ],
        "verification_boundary": "Exact determinant weights, exhaustive proper cycle inventory, overlapping forest identity, radial moments and local-polynomial annihilation; independent graph, matrix, momentum-region and Cauchy tests. Analytic continuation and convergence arguments are written proofs, not formalization or independent peer review. Native, ordinary, CLI and direct science retain unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The quadratic-forest report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.145.PAIRED_QUADRATIC_NONLOCAL_REMAINDERS replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
