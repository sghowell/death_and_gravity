"""Read-only complete matched GY14 two-loop Phi forward-coefficient certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_full_phi_normalization import verify as parent

from . import audit, bounds, calibration, contact, matching, ownership

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-full-two-loop-amplitude.json"
PARENT_SHA = "a75a1b9dc9f16e98e9d58500b3b36cb48f17d2744f3973c95108234facf37080"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_full_two_loop_amplitude/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen complete two-loop amplitude parent report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_156_fully_rebuilt": PARENT_SHA,
        "same_GY14_nongravitational_reference_model": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
        "complete_fixed_order_Phi_pole_and_b2_not_vacuum_truncation_or_V_G_B": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A complete two-loop amplitude proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.157.COMPLETE_MATCHED_GY14_TWO_LOOP_PHI_FORWARD_COEFFICIENT",
        "date": "2026-09-10",
        "status": "COMPLETE_MATCHED_FIXED_ORDER_GY14_TWO_LOOP_PHI_B2_AND_UNIT_DISC_POLE; NOT_VACUUM_SOURCE_PHYSICAL_TRUNCATION_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/ownership.md",
            "notes/matching.md",
            "notes/contact.md",
            "notes/amplitude.md",
            "notes/bounds.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_four_point_ownership_and_coordinate_matching": serialize(
            {
                "ownership": payload(ownership.data()),
                "matching": payload(matching.data()),
                "contact": payload(contact.data()),
            }
        ),
        "formal_two_loop_assembly_and_error_rule": serialize(payload(bounds.data())),
        "actual_complete_canonical_Phi_amplitude_enclosure": serialize(
            payload(calibration.data())
        ),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "All 192 scalar interaction refinements and all four fermionic quartic families, including their assigned proper and overall MS references, are assembled in the physical-Phi hybrid representation. The full regulated G/L/M field map, induced M2 and fundamental G1 square are retained. The first coordinate variation of the finite hybrid one-loop amplitude has scalar weight four and fermion weight two. Pointwise contact cancellation is composed without assuming commuting parameter directions. The complete two-loop relative b2 bound is below 1e-7; the complete one-plus-two-loop bound is below 1e-6, so the formal coefficient is positive. The complete unit-disc canonical pole is inherited unchanged. Vacuum/source references, physical truncation, V contours/cuts, G, B and original P8 remain open.",
        "not_established": [
            "The complete vacuum/source reference assembly or global quantum potential",
            "A physical higher-loop or finite-EFT truncation error",
            "External-fermion second-order canonical dictionaries beyond this Phi result",
            "V contour/cut control, finite-gravity G, common-parent B or original P8 closure",
        ],
        "verification_boundary": "Complete explicit four-point family and counterterm ownership, regulator-complete inherited bare maps, new finite-coordinate homogeneity and exact rational assembly. Independent tests compose full mapped amplitudes, extract the forward coefficient on a complex circle, retain induced heavy-mass and cubic-square terms, and test noncommuting-contact cancellation. Written analytic arguments and exact replay are not formalization or independent peer review. Native, ordinary, CLI and direct science retain unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete two-loop amplitude report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.157.COMPLETE_MATCHED_GY14_TWO_LOOP_PHI_FORWARD_COEFFICIENT replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
