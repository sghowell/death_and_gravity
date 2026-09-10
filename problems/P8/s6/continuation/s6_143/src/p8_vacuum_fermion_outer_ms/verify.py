"""Read-only MS fermion outer-ms certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_fermion_vertex_chord import verify as parent

from . import audit, calibration, conversion, dimensional, laurent, remainder

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-fermion-outer-ms.json"
PARENT_SHA = "2b0f5a84ac0e2fc1a0a5ffaba7c4930c4a0006c9d5d7cc66853475b21db33ea0"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_fermion_outer_ms/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen MS fermion outer-ms report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_142_fully_rebuilt": PARENT_SHA,
        "same_GY14_nongravitational_reference_model": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
        "only_one_paired_family_outer_reference_advances": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A MS fermion outer-ms gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.143.FERMION_OUTER_MS_REFERENCE_CONVERSION",
        "date": "2026-09-10",
        "status": "FULL_HEAVY_VERTEX_FAMILY_WITH_OUTER_MS_REFERENCE_CONVERSION; NOT_FULL_MATCHED_TWO_LOOP_OTHER_CONVERSIONS_HIGHER_LOOP_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/dimensional.md",
            "notes/finite_part.md",
            "notes/conversion.md",
            "notes/remainder.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "finite_outer_MS_reference": serialize(
            {"laurent": payload(laurent.data()), "remainder": payload(remainder.data())}
        ),
        "dimensional_OS_spectral_representation": serialize(
            payload(dimensional.data())
        ),
        "local_parent_reference_conversion": serialize(payload(conversion.data())),
        "actual_converted_family_frontier": serialize(
            {
                "actual": payload(calibration.data()),
                "primitive_frontier": audit.frontier(),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete dimensional on-shell fermion insertion defines the outer zero-momentum reference before MS subtraction. Its leading gamma-function integral has finite coefficient 47/18+pi^2/12. The exact nonzero mass-ratio correction retains its simple pole and finite prefactor products; its finite part is bounded by r(20-2log r). A literal local-parent parameter map fixes the conversion sign and b2 coefficient. The S6.138 full heavy-vertex family retains a b2 bound below 1e-622 in the common MS interaction scheme. Other primitive rows, finite field/parameter conversions, complete two-loop error, V/G/B and original P8 remain open.",
        "not_established": [
            "Five other wholly unevaluated primitive fermion-sector rows",
            "Other finite field, parameter and counterterm matching contributions",
            "Complete two-loop amplitude, canonical pole or matched error",
            "Higher-loop truncation control, V contours, finite-gravity G, common-parent B or original P8 closure",
        ],
        "verification_boundary": "Exact dimensional OS and spectral normalization, Gamma/polygamma and Laurent identities, an integrable finite remainder, literal parent-reference conversion and mutation-checked scope. Independent tests include the nonzero-regulator radial integral and complex-regulator finite-part extraction. Analytic continuation and regulator-limit arguments are written proofs, not proof-assistant formalization or independent peer review. Native, ordinary, CLI and direct science use unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The MS fermion outer-ms differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.143.FERMION_OUTER_MS_REFERENCE_CONVERSION replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
