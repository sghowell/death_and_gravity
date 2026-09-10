"""Read-only once-fixed finite potential contact and its two-loop insertion."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_two_loop_wineglass import verify as parent

from . import audit, calibration, insertion, potential, radial, subtraction

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-two-loop-finite-contact.json"
PARENT_SHA = "148c12cdb2fc295d8cf7ac3fabfb2e3a5d21807de25bbbcf8a6ef0bf15412e37"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_two_loop_finite_contact/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen subtracted wineglass report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_124_fully_rebuilt": PARENT_SHA,
        "same_polynomial_Hessian_actual_couplings_and_quartic_potential_condition": True,
        "same_full_I0_reference_and_complete_one_loop_amplitude": True,
        "all_four_disjoint_raw_graph_groups_retained_without_recounting": True,
        "no_full_normalization_all_orders_or_original_P8_closure_inferred": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = {k: bool(v) for k, v in audit.gates().items()}
    if not all(v is True for v in gates.values()):
        raise ValueError("A fixed finite-contact insertion gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.125.POLYNOMIAL_VACUUM_TWO_LOOP_FIXED_FINITE_CONTACT",
        "date": "2026-09-10",
        "status": "COMPLETE_ONCE_FIXED_FINITE_POTENTIAL_QUARTIC_CONTACT_FULL_HESSIAN_ANCHORED_INTEGRALS_AND_TWO_LOOP_INSERTION_WITH_LINEAR_INHERITED_LOCAL_REFERENCE; NOT_COMPLETE_TWO_LOOP_POLE_LSZ_SOURCE_MATCHING_ALL_HIGHER_LOOP_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/potential.md",
            "notes/radial.md",
            "notes/insertion.md",
            "notes/bounds.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": len(rows),
        "proof_checks": gates,
        "full_Hessian_and_fixed_quartic_contact_condition": serialize(
            payload(potential.data())
        ),
        "independent_anchored_finite_contact_integrals": serialize(
            payload(radial.data())
        ),
        "complete_insertion_and_same_local_reference": serialize(
            {
                "insertion": payload(insertion.data()),
                "subtraction": payload(subtraction.data()),
            }
        ),
        "actual_fixed_contact_and_known_error_budget": serialize(
            {
                "actual": payload(calibration.data()),
                "zero_radial_value": calibration.point(0),
                "unit_radial_value": calibration.point(1),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The original constant-field quartic condition fixes a strictly negative finite contact after the full I0 reference subtraction. Its exact value follows from the literal two-field Hessian and two independent anchored heavy-mass integrals. Its two-loop insertion is the fixed-g,M quartic variation of the complete renormalized one-loop amplitude, including both heavy triangles and the linear inherited local reference counterterms. The insertion has absolute forward-coefficient upper bound below 7 times 10^-612 and below 2 times 10^-12 of tree. Adding all four raw-graph groups leaves the known two-loop contributions below 10^-7 of tree; the known one-plus-two-loop contributions remain below 10^-6 of tree. Full quantitative pole/residue/LSZ and source-aware matching are not inferred from these contribution bounds.",
        "not_established": [
            "The signed or exact value of the complete two-loop forward coefficient",
            "Complete quantitative two-loop light pole residue LSZ and source-aware matching",
            "All higher-loop truncation errors all-energy contours or a nonperturbative ultraviolet completion",
            "Bare ultraviolet stability inferred from a small finite counterterm",
            "Finite-gravity IR/Regge Delta absolute coupled renormalization or a common controlled bounce parent",
            "Full V/G/B classification any excluded original P8 row or original P8 closure",
        ],
        "verification_boundary": "Native literal full Hessian and frozen-kernel identity, fixed potential Taylor cancellation, independently anchored positive heavy-mass integrals, complete fixed-parameter quartic variation, local reference counterterms and exact rational actual bounds. Continuous sign, integration, renormalization and Cauchy arguments are source-pinned written proofs, not proof-assistant formalized, nonperturbative or peer reviewed. Every ancestor, own source and report field is rebuilt read-only with unmodified scientific SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The fixed finite-contact report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.125.POLYNOMIAL_VACUUM_TWO_LOOP_FIXED_FINITE_CONTACT replay passed; full V/G/B and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
