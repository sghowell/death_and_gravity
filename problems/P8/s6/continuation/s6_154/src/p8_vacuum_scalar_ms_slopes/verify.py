"""Read-only complete scalar quadratic finite-MS-slope certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_scalar_insertion_ms import verify as parent

from . import audit, bounds, calibration, cauchy, conversion, sunset

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-scalar-ms-slopes.json"
PARENT_SHA = "46311af32da0d789a75132b6c9515c9e34f04966a4422deb794e31bd35e34a0d"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_scalar_ms_slopes/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen scalar MS slope parent report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_153_fully_rebuilt": PARENT_SHA,
        "same_GY14_nongravitational_reference_model": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
        "only_scalar_quadratic_finite_MS_reference_and_converted_OS_bound_advance": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A scalar quadratic MS proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.154.COMPLETE_SCALAR_QUADRATIC_MS_SLOPE_REFERENCE",
        "date": "2026-09-10",
        "status": "COMPLETE_SCALAR_QUADRATIC_FINITE_MS_SLOPE_AND_CONVERTED_OS_BOUND; NOT_FULL_GY14_FIELD_COUPLING_MATCHING_TRUNCATION_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/sunset.md",
            "notes/forests.md",
            "notes/slopes.md",
            "notes/pole.md",
            "notes/reference.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "dimensional_local_sunset_and_complete_Cauchy_bounds": serialize(
            {"sunset": payload(sunset.data()), "cauchy": payload(cauchy.data())}
        ),
        "proper_MS_conversion_and_integral_enclosures": serialize(
            {"conversion": payload(conversion.data()), "bounds": payload(bounds.data())}
        ),
        "actual_finite_MS_slope_and_converted_OS_enclosure": serialize(
            payload(calibration.data())
        ),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The finite scalar two-loop MS slope is bounded for all 32 quadratic refinements and required counterterm insertions, retaining the fixed inner physical OS convention and interaction-MS references. The local sunset slope has pole L^2/(24Q^2 epsilon); a convergent simplex representation fixes its finite part before outer OS. First-sheet integrated Cauchy bounds control the remaining slopes, including the nested alpha_0^2 term. The proper-reference conversion and squared-heavy geometric remainder are retained. The complete scalar finite MS slope is below 1e-18, and the converted one-plus-two-loop OS inverse retains a unique unit-residue pole on the unit disc after imposing the outer physical conditions. This intermediate scalar statement is not a full GY14 finite field/coupling map, physical higher-loop error, global pole theorem or V/G/B closure.",
        "not_established": [
            "Global finite field/coupling map, evanescent re-expansion and cross terms",
            "Complete GY14 vacuum/source and canonical amplitude/pole assembly",
            "Physical higher-loop truncation, V contours/cuts or finite-gravity G",
            "Common-parent B or original P8 closure",
        ],
        "verification_boundary": "Exact Schwinger/projective identities, regulated derivative and proper-reference algebra, the full inherited graph inventory, proven first-sheet integrated majorants and positive parameter bounds. Independent tests probe nonzero regulators, complex finite-part extraction, direct kernel derivatives and geometric remainders, with input and scope mutations. Written analytic arguments and exact replay are not formalization or independent peer review. Native, ordinary, CLI and direct science retain unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The scalar quadratic MS report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.154.COMPLETE_SCALAR_QUADRATIC_MS_SLOPE_REFERENCE replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
