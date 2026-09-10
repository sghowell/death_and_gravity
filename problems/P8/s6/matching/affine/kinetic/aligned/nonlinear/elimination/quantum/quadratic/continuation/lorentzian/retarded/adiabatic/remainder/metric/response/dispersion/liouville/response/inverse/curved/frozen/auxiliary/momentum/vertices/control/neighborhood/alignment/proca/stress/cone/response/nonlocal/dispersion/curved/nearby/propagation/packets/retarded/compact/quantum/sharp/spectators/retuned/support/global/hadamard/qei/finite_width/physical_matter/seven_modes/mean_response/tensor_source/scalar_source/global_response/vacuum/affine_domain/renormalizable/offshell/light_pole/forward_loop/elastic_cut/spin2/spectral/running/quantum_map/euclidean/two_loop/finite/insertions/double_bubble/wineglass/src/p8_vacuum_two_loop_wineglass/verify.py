"""Read-only complete wineglass group and combined raw-refinement bound."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_two_loop_double_bubble import verify as parent

from . import audit, calibration, integrals, kernel, selection, subtraction

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-two-loop-wineglass.json"
PARENT_SHA = "96e5cd10cd9d7749dfdd71264e4edcfedc0eb690e291cc94fd2f423c0593466b"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_two_loop_wineglass/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen subtracted double-bubble report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_123_fully_rebuilt": PARENT_SHA,
        "same_actual_graphs_couplings_and_first_sheet": True,
        "same_proper_I0_and_explicit_recursive_local_reference_framework": True,
        "previous_three_disjoint_raw_graph_groups_retained": True,
        "complete_raw_graph_count_not_full_two_loop_amplitude_or_P8_closure": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = {k: bool(v) for k, v in audit.gates().items()}
    if not all(v is True for v in gates.values()):
        raise ValueError("A subtracted wineglass gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.124.POLYNOMIAL_VACUUM_TWO_LOOP_SUBTRACTED_WINEGLASS",
        "date": "2026-09-10",
        "status": "COMPLETE_SIXTEEN_ACTUAL_WINEGLASS_REFINEMENTS_RECURSIVE_LOCAL_SUBTRACTIONS_TWELVE_COMPLEX_ROUTINGS_DECAYING_LOGARITHM_DIFFERENCE_AND_ALL_RADIUS_BOUND_ALL_192_RAW_REFINEMENTS_REPRESENTED; NOT_FINITE_POTENTIAL_INSERTIONS_COMPLETE_TWO_LOOP_POLE_LSZ_MATCHING_HIGHER_LOOP_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/selection.md",
            "notes/subtraction.md",
            "notes/kernel.md",
            "notes/bounds.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": len(rows),
        "proof_checks": gates,
        "exact_remaining_refinements_and_channel_weights": serialize(
            payload(selection.data())
        ),
        "actual_restricted_forests_and_recursive_local_reference": serialize(
            payload(subtraction.data())
        ),
        "all_routings_decaying_logarithm_and_integrated_majorants": serialize(
            {
                "kernel": payload(kernel.data()),
                "integrals": payload(integrals.data()),
            }
        ),
        "actual_wineglass_and_complete_raw_group_bound": serialize(
            {
                "actual": payload(calibration.data()),
                "zero_radial_enclosure": calibration.point(0),
                "unit_radial_enclosure": calibration.point(1),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The remaining 16 raw wineglass refinements have their actual proper and overall restricted forests evaluated with the fixed full I0 reference and an explicit recursive zero-momentum extension. All 12 external-label routings keep the inner logarithm on its first sheet. An anchored logarithm difference decays at large radius, while every heavy-weighted remainder retains the complete propagators. Exact radial primitives and actual rational comparisons bound this family's second forward coefficient below 10^-611 and below 3 times 10^-12 of tree. The four disjoint graph groups represent all 192 bare refinements with their specified subtractions and remain below 10^-7 of tree. The signed correction, fixed finite-potential insertions and complete two-loop normalization are not determined by this raw-graph count.",
        "not_established": [
            "The signed or exact value of this subtracted family or the combined graph groups",
            "Fixed lower-order finite-potential counterterms inserted into one-loop amplitudes",
            "Complete two-loop local counterterm light mass residue LSZ and source-aware matching",
            "All higher-loop truncation errors all-energy contours or an ultraviolet completion",
            "Finite-gravity IR/Regge Delta absolute coupled renormalization or a common controlled bounce parent",
            "Full V/G/B classification any excluded original P8 row or original P8 closure",
        ],
        "verification_boundary": "Native complete graph and external-assignment selection, actual restricted-forest terms and local parameter variations, all twelve bilinear routings, anchored logarithmic and rational primitives, convergence majorants and exact rational actual bounds. Continuous contour continuation, common-regulator subtraction, absolute integration and Cauchy estimates are source-pinned written proofs, not proof-assistant formalized, nonperturbative or peer reviewed. Every ancestor, own source and report field is rebuilt read-only with unmodified scientific SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The subtracted wineglass report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.124.POLYNOMIAL_VACUUM_TWO_LOOP_SUBTRACTED_WINEGLASS replay passed; full V/G/B and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
