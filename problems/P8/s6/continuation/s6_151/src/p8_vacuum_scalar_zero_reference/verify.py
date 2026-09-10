"""Read-only equal-mass reference and wineglass MS conversion certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_finite_contact_conversion import verify as parent

from . import anchors, audit, bounds, calibration, conversion, sector

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-scalar-zero-reference.json"
PARENT_SHA = "4f3154ca66df6c4fd39426bd87081f1f15078e687822908c2e718f80ca00e673"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_scalar_zero_reference/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen scalar zero-reference parent report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_150_fully_rebuilt": PARENT_SHA,
        "same_GY14_nongravitational_reference_model": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
        "only_wineglass_interaction_forest_conversion_advances": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A wineglass MS proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.151.COMPLETE_WINEGLASS_MS_INTERACTION_FOREST",
        "date": "2026-09-10",
        "status": "COMPLETE_WINEGLASS_INTERACTION_FOREST_MS_WITH_REGULATED_LOCAL_ANCHOR_AND_FULL_SCALE_TERM; NOT_FULL_SCALAR_GY14_MATCHING_CANONICAL_POLE_TRUNCATION_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/sector.md",
            "notes/anchors.md",
            "notes/conversion.md",
            "notes/bounds.md",
            "notes/reference.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "regulated_sector_and_local_MS_anchors": serialize(
            {"sector": payload(sector.data()), "anchors": payload(anchors.data())}
        ),
        "complete_wineglass_conversion_and_bounds": serialize(
            {"conversion": payload(conversion.data()), "bounds": payload(bounds.data())}
        ),
        "actual_wineglass_MS_enclosure": serialize(payload(calibration.data())),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete sixteen-refinement wineglass interaction forest is converted to MS at mu=mF at fixed canonical reference coordinates. A six-sector regulator-first equal-mass vacuum derivative fixes both the old recursive zero reference and the distinct proper-MS anchor, retaining finite epsilon-times-pole terms. The actual conversion also includes the full finite inner-scale change over all outer light/heavy propagators. Its local reference has the exact heavy-tree b2 variation. The complete family is bounded below 3e-611 and below 8e-12 relative to tree. Primitive statuses are unchanged. Other scalar families, full parameter/field re-expansion and cross terms, vacuum/canonical assembly, physical truncation and V/G/B remain open.",
        "not_established": [
            "Other scalar interaction-family MS conversions and complete finite field/coupling map",
            "Complete GY14 vacuum/source references, matched amplitude or canonical pole/error",
            "Physical higher-loop truncation, V contour/cut control or finite-gravity G",
            "Common-parent B or original P8 closure",
        ],
        "verification_boundary": "Exact six-sector endpoint subtraction, equal-mass vacuum differentiation and Laurent products, a convergent bounded finite constant, complete regulated forest algebra and all-radius full-propagator bounds. Independent tests use the bubble/radial hypergeometric representation, nonzero regulators, finite-part extraction and mutation checks. Written analytic arguments and exact replay are not formalization or independent peer review. Native, ordinary, CLI and direct science retain unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The wineglass MS report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.151.COMPLETE_WINEGLASS_MS_INTERACTION_FOREST replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
