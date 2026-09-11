"""Read-only complete two-loop GY14 vacuum-reference certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_full_two_loop_amplitude import verify as parent

from . import audit, bounds, calibration, forest, sector, source

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-full-vacuum-reference.json"
PARENT_SHA = "ca620649aa67631a60681b9e32f2cac204cae0f530230261a9490d1064b5651f"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_full_vacuum_reference/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen complete vacuum-reference parent report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_157_fully_rebuilt": PARENT_SHA,
        "same_GY14_nongravitational_reference_model": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
        "complete_vacuum_reference_not_full_second_source_truncation_or_V_G_B": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A complete vacuum-reference proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.158.COMPLETE_MATCHED_TWO_LOOP_VACUUM_REFERENCE_AND_SOURCE_SQUARE",
        "date": "2026-09-10",
        "status": "COMPLETE_MATCHED_GY14_FIXED_ORDER_TWO_LOOP_VACUUM_REFERENCE_AND_REGULATED_FIRST_SOURCE_SQUARE; NOT_FULL_SECOND_SOURCE_GLOBAL_POTENTIAL_PHYSICAL_TRUNCATION_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/forest.md",
            "notes/sector.md",
            "notes/source.md",
            "notes/reference.md",
            "notes/bounds.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_vacuum_forest_and_regulated_source_ownership": serialize(
            {"forest": payload(forest.data()), "source": payload(source.data())}
        ),
        "sunset_analytic_continuation_and_finite_part_bounds": serialize(
            {"sector": payload(sector.data()), "bounds": payload(bounds.data())}
        ),
        "actual_complete_vacuum_reference_enclosure": serialize(
            payload(calibration.data())
        ),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The exact Gaussian-heavy source square, its induced mass insertion and the heavy-reducible vacuum graph cancel only as a common-regulator triple. The remaining complete scalar vacuum forest includes physical Phi mass/residue and heavy-MS-mass counterterms. Six corner-subtracted sectors and a complex-regulator Cauchy estimate bound its finite part. Adding both already assigned fermionic vacuum forests, with all 42 gauge-coupled states, gives the complete two-loop vacuum-reference bound below 1e595 and below 1e-203 of the complete first-loop vacuum lower bound. The physical vacuum-energy-zero condition fixes the corresponding first and second finite local references. The first Phi field map leaves the regulated first-loop determinant invariant. The full second H-source, global potential, physical truncation and V/G/B remain open.",
        "not_established": [
            "The complete second H-source reference, although its square does not enter this order",
            "A global effective potential, nonperturbative gauge vacuum or semiclassical backreaction theorem",
            "A physical finite-EFT higher-order truncation bound",
            "V contours/cuts, finite-gravity G, common-parent B or original P8 closure",
        ],
        "verification_boundary": "Explicit Wick and proper-counterterm ownership, exact regulated Gaussian elimination, common-regulator finite products, six-sector analytic continuation, Cauchy coefficient bounds and exact rational actual-parameter power caps. Independent tests check literal Gaussian Wick pairings, source/pole products, simplex changes of variables, zero-dimensional sunset normalization, regulator-circle estimates and scope controls. Written analytic arguments and exact replay are not formalization or independent peer review. Native, ordinary, CLI and direct science retain unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete vacuum-reference report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.158.COMPLETE_MATCHED_TWO_LOOP_VACUUM_REFERENCE_AND_SOURCE_SQUARE replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
