"""Read-only local forest subtraction and actual double-bubble group bound."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_two_loop_insertions import verify as parent

from . import audit, bubble, calibration, forests, selection, triangle

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-two-loop-double-bubble.json"
PARENT_SHA = "dbe31886cc72175d68767af86c742d9fa526087241039034e081f21db0c2146c"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_two_loop_double_bubble/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen grouped insertion report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_122_fully_rebuilt": PARENT_SHA,
        "same_actual_couplings_and_complete_heavy_graph_inventory": True,
        "same_regulated_I0_reference_and_first_sheet_routing": True,
        "previous_two_disjoint_groups_retained": True,
        "no_complete_two_loop_or_original_P8_result_inferred": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = {k: bool(v) for k, v in audit.gates().items()}
    if not all(v is True for v in gates.values()):
        raise ValueError("A subtracted double-bubble gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.123.POLYNOMIAL_VACUUM_TWO_LOOP_SUBTRACTED_DOUBLE_BUBBLE",
        "date": "2026-09-09",
        "status": "COMPLETE_TWENTY_FOUR_ACTUAL_DOUBLE_BUBBLE_REFINEMENTS_RECURSIVE_LOCAL_ZERO_MOMENTUM_FORESTS_FIXED_REFERENCE_FINITE_FACTOR_PRODUCTS_AND_INTEGRATED_COMPACT_FORWARD_BOUND; NOT_REMAINING_SIXTEEN_WINEGLASS_FINITE_POTENTIAL_INSERTIONS_COMPLETE_TWO_LOOP_POLE_LSZ_HIGHER_LOOP_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/selection.md",
            "notes/forests.md",
            "notes/factors.md",
            "notes/bounds.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": len(rows),
        "proof_checks": gates,
        "exact_disjoint_refinement_selection": serialize(payload(selection.data())),
        "every_actual_restricted_forest_and_local_subtraction": serialize(
            payload(forests.data())
        ),
        "complete_finite_bubble_and_heavy_triangle_factors": serialize(
            {
                "bubble": payload(bubble.data()),
                "triangle": payload(triangle.data()),
            }
        ),
        "actual_integrated_group_and_combined_partial_bound": serialize(
            {
                "actual": payload(calibration.data()),
                "zero_channel_diagnostic": calibration.point(0),
                "forward_center_diagnostic": calibration.point(2),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "Every actual restricted forest of the 24 subtraction-dependent factorizing double-bubble refinements is evaluated in the same regulated one-loop reference with an explicit recursive zero-momentum extension for shared two-loop cores. The shared overall subtraction, not an illegal pair of overlapping cores, supplies the missing reference square. Local full-model parameter variations and permitted iterated propagator counterterms implement the result. The selected renormalized group factorizes into finite light bubbles and complete heavy triangles, excluding eight already bounded finite refinements. Its integrated forward coefficient has absolute upper bound below 2 times 10^-614 and below 3 times 10^-15 of tree. Three disjoint groups represent 176 raw refinements and remain below 10^-7 of tree, without determining the signed correction.",
        "not_established": [
            "The signed exact contribution of this subtracted graph family",
            "The other 16 bare wineglass refinements and their integrated subtractions",
            "Lower-order finite-potential counterterm insertions and complete two-loop counterterm accounting",
            "Complete two-loop mass residue LSZ or source-aware quantum field-map control",
            "All higher-loop errors all-energy contours finite-gravity IR/Regge Delta or a common controlled bounce parent",
            "Full V/G/B classification any excluded original P8 row or original P8 closure",
        ],
        "verification_boundary": "Native exhaustive graph selection and actual restricted-forest sums, local full-model parameter identities, anchored finite logarithm and all-radius triangle integrals, exact symmetry factors and rational actual bounds. Continuous analytic continuation, regulator removal, recursive local renormalization and Cauchy estimates are source-pinned written proofs, not proof-assistant formalized, nonperturbative or peer reviewed. Every ancestor, own source and report field is rebuilt read-only with unmodified scientific SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The subtracted double-bubble report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.123.POLYNOMIAL_VACUUM_TWO_LOOP_SUBTRACTED_DOUBLE_BUBBLE replay passed; full V/G/B and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
