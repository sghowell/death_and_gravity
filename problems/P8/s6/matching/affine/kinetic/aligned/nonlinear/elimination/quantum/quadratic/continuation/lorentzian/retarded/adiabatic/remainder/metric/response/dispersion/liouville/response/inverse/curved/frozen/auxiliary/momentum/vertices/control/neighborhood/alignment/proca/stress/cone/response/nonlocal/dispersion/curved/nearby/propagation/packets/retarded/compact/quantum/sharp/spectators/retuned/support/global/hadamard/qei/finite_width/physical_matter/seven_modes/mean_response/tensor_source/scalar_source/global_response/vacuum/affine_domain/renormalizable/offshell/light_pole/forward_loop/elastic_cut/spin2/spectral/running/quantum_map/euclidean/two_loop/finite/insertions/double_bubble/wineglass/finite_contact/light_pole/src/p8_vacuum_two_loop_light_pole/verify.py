"""Read-only complete canonical two-loop light-pole checkpoint."""

import argparse
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_two_loop_finite_contact import verify as parent

from . import analytic, audit, calibration, denominators, graphs, sectors, subtraction

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-two-loop-light-pole.json"
PARENT_SHA = "80ea67b9827236b172bc7f639988d293f64eac9a7b27de0d80533f7afcb4ace8"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_two_loop_light_pole/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen fixed finite-contact report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_125_fully_rebuilt": PARENT_SHA,
        "same_actual_polynomial_vacuum_and_canonical_parameters": True,
        "same_once_fixed_potential_contact_and_entire_I0_reference": True,
        "same_inner_OS_light_kernel_and_canonical_interaction_convention": True,
        "no_complete_four_point_source_map_or_original_P8_closure_inferred": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = {k: bool(v) for k, v in audit.gates().items()}
    if not all(v is True for v in gates.values()):
        raise ValueError("A quantitative two-loop light-pole gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.126.POLYNOMIAL_VACUUM_TWO_LOOP_LIGHT_POLE",
        "date": "2026-09-10",
        "status": "COMPLETE_CANONICAL_TWO_LOOP_LIGHT_TWO_POINT_INVENTORY_LOCAL_SUBTRACTIONS_INTEGRATED_FIRST_SHEET_BOUND_AND_UNIT_DISC_POLE_RESIDUE; NOT_COMPLETE_TWO_LOOP_FOUR_POINT_SOURCE_MATCHING_ALL_HIGHER_LOOP_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/inventory.md",
            "notes/denominators.md",
            "notes/sectors.md",
            "notes/subtraction.md",
            "notes/analytic.md",
            "notes/bounds.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": len(rows),
        "proof_checks": gates,
        "complete_two_point_graph_and_forest_inventory": serialize(
            payload(graphs.data())
        ),
        "independent_denominators_and_integrated_sectors": serialize(
            {
                "denominators": payload(denominators.data()),
                "sectors": payload(sectors.data()),
            }
        ),
        "actual_subtractions_and_complex_analytic_bounds": serialize(
            {
                "subtractions": payload(subtraction.data()),
                "analytic": payload(analytic.data()),
            }
        ),
        "actual_two_loop_light_pole_and_error_budget": serialize(
            {
                "actual": payload(calibration.data()),
                "unit_disc": calibration.point(1),
                "half_disc": calibration.point(Fraction(1, 2)),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "All 32 reduced Gaussian-heavy two-loop light two-point refinements and the required inherited one-loop local insertions are grouped in the fixed canonical scheme. Independent spanning-tree/two-forest and Laplacian/adjugate constructions give explicit first-sheet denominator squares. All finite sunset sectors, the twice-differentiated local sunset and the subtracted mixed, same-heavy and nested contributions have complete integrated bounds. The actual two-loop OS self-energy obeys |Pi2_R(s)| < 10^-18 |s-1|^2 on |s-1|<=1; adding the inherited one-loop bound preserves this strict bound. The inverse propagator through two loops therefore has exactly the mass-one pole in that disc, with unit residue. This is a fixed-order canonical two-point result, not a nonperturbative propagator, complete two-loop four-point/source-matching calculation or original P8 closure.",
        "not_established": [
            "The exact or signed two-loop self-energy, a global pole count or nonperturbative spectral positivity",
            "An independent complete canonical two-loop four-point normalization ledger",
            "Two-loop source-aware derivative matching or controlled cosmological loops",
            "All higher-loop truncation errors all-energy contours or a UV completion",
            "Finite-gravity IR/Regge Delta absolute coupled renormalization or a common controlled bounce parent",
            "Full V/G/B classification any excluded original P8 row or original P8 closure",
        ],
        "verification_boundary": "Native exhaustive graph/UV-forest inventory, independent denominator polynomials and current sum-of-squares identities, all ordered sectors and anchored elementary integrals, actual inherited counterterm cancellations and exact rational parameter bounds. Continuous first-sheet continuation, logarithm and integration bounds, on-shell subtraction and pole arguments are source-pinned written proofs, not proof-assistant formalized, nonperturbative or peer reviewed. Every ancestor, own source and report field is rebuilt read-only with unmodified scientific SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The two-loop light-pole report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.126.POLYNOMIAL_VACUUM_TWO_LOOP_LIGHT_POLE replay passed; full V/G/B and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
