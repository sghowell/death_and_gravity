"""Read-only complete double-bubble interaction-forest MS certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_scalar_zero_reference import verify as parent

from . import audit, bounds, calibration, conversion, forests, forward

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-double-bubble-ms.json"
PARENT_SHA = "347356124e5cda2c5d2c68510c5c2f7325414f0559b49219c068b5d8e81295ad"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_double_bubble_ms/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen double-bubble MS parent report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_151_fully_rebuilt": PARENT_SHA,
        "same_GY14_nongravitational_reference_model": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
        "only_double_bubble_interaction_forest_conversion_advances": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A double-bubble MS proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.152.COMPLETE_DOUBLE_BUBBLE_MS_INTERACTION_FOREST",
        "date": "2026-09-10",
        "status": "COMPLETE_DOUBLE_BUBBLE_INTERACTION_FOREST_MS_WITH_ACTUAL_NESTED_FORESTS_AND_BOTH_SCALE_TERMS; NOT_FULL_SCALAR_GY14_MATCHING_CANONICAL_POLE_TRUNCATION_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/forests.md",
            "notes/factors.md",
            "notes/conversion.md",
            "notes/bounds.md",
            "notes/reference.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "actual_MS_forests_and_full_scale_factors": serialize(
            {
                "forests": payload(forests.data()),
                "conversion": payload(conversion.data()),
            }
        ),
        "exact_quadratic_forward_term_and_bounds": serialize(
            {"forward": payload(forward.data()), "bounds": payload(bounds.data())}
        ),
        "actual_double_bubble_MS_enclosure": serialize(payload(calibration.data())),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete 24-refinement double-bubble interaction forest is converted to MS at mu=mF at fixed canonical reference coordinates. Every actual single, disjoint and overlapping forest is evaluated before regulator removal. Raw and nested overall projections retain their distinct finite-bubble pole terms, and the forbidden shared pair remains absent. The resulting full scale shift includes both heavy triangles and both linear and quadratic terms. The latter retains its exact heavy-tree forward coefficient, including repeated ordinary heavy insertions. The complete assigned family is below 7e-612 and below 2e-12 relative to tree. Other scalar families, full field/coupling re-expansion and cross terms, vacuum/canonical assembly, physical truncation and V/G/B remain open.",
        "not_established": [
            "Other scalar interaction-family MS conversions and complete finite field/coupling map",
            "Complete GY14 vacuum/source references, matched amplitude or canonical pole/error",
            "Physical higher-loop truncation, V contour/cut control or finite-gravity G",
            "Common-parent B or original P8 closure",
        ],
        "verification_boundary": "Exact actual graph-forest projections, regulator-first cancellation, complete full-propagator factorization, exact heavy-tree forward coefficients and all-radius bounds. Independent tests use explicit Laurent projections, all actual forest classes, dimensional bubbles, complex-disc integrals, contour coefficients and mutation controls. Written analytic arguments and exact replay are not formalization or independent peer review. Native, ordinary, CLI and direct science retain unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The double-bubble MS report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.152.COMPLETE_DOUBLE_BUBBLE_MS_INTERACTION_FOREST replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
