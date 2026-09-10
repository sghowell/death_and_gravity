"""Read-only complete scalar OS insertion interaction-forest MS certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_double_bubble_ms import verify as parent

from . import audit, bounds, calibration, forward, inner, outer

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-scalar-insertion-ms.json"
PARENT_SHA = "df058b81fec53d4d6d2ee20c83935a59694bda97dcbabece17ad0ab1fc785505"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_scalar_insertion_ms/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen scalar insertion MS parent report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_152_fully_rebuilt": PARENT_SHA,
        "same_GY14_nongravitational_reference_model": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
        "only_scalar_OS_insertion_conversion_advances_and_raw_families_assembled": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A scalar insertion MS proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.153.COMPLETE_SCALAR_INSERTION_MS_INTERACTION_FOREST",
        "date": "2026-09-10",
        "status": "COMPLETE_SCALAR_OS_INSERTION_INTERACTION_FOREST_MS_AND_FOUR_RAW_FAMILY_BOUND; NOT_FULL_SCALAR_GY14_MATCHING_CANONICAL_POLE_TRUNCATION_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/inner.md",
            "notes/outer.md",
            "notes/conversion.md",
            "notes/bounds.md",
            "notes/reference.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "dimensional_inner_slope_and_outer_MS_reference": serialize(
            {
                "inner": payload(inner.data()),
                "outer": payload(outer.data()),
            }
        ),
        "exact_forward_reference_and_bounds": serialize(
            {"forward": payload(forward.data()), "bounds": payload(bounds.data())}
        ),
        "actual_insertion_and_four_scalar_family_enclosure": serialize(
            payload(calibration.data())
        ),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete 64-refinement scalar OS insertion interaction forest is converted to MS at mu=mF at fixed canonical reference coordinates. The full dimensional slope fixes the finite outer reference F_alpha=(g/Q^2) integral b(2ell-logD), retaining the epsilon slope coefficient in its outer pole product. The inner physical mass/residue grouping and both convergent decaying insertions are unchanged. The exact local heavy-tree variation bounds the conversion below 1e-1008 and the complete insertion family below 1e-614. Together with the unchanged 88 UV-finite refinements and the preceding double-bubble and wineglass conversions, all four raw scalar interaction families remain below 3e-607 and below 1e-7 relative to tree. This is not global MS field/coupling matching, the full GY14 canonical amplitude/pole, vacuum/source assembly, physical truncation or V/G/B.",
        "not_established": [
            "Global finite field/coupling map and its cross terms, beyond the four raw scalar interaction families",
            "Complete GY14 vacuum/source references, matched amplitude or canonical pole/error",
            "Physical higher-loop truncation, V contour/cut control or finite-gravity G",
            "Common-parent B or original P8 closure",
        ],
        "verification_boundary": "Exact dimensional on-shell differentiation, regulator-first slope/pole products, complete two-position insertion grouping, exact three-channel heavy-tree variation and positive parameter-integral bounds. Independent tests use nonzero-regulator derivatives, complex finite-part extraction, large-momentum asymptotics, contour coefficients and mutation controls. Written analytic arguments and exact replay are not formalization or independent peer review. Native, ordinary, CLI and direct science retain unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The scalar insertion MS report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.153.COMPLETE_SCALAR_INSERTION_MS_INTERACTION_FOREST replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
