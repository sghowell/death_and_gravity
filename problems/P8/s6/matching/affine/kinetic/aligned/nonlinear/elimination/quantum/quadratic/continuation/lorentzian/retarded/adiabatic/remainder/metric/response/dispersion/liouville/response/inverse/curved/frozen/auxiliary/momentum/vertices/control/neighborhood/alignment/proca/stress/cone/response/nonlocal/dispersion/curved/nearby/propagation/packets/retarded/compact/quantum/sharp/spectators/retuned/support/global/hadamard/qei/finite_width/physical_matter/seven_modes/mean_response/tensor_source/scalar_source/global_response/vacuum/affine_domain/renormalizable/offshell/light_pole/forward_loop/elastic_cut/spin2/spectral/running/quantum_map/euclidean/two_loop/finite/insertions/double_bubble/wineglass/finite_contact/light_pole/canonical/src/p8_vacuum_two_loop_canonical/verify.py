"""Read-only complete canonical two-loop forward-amplitude checkpoint."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_two_loop_light_pole import verify as parent

from . import assembly, audit, calibration, contractions, fields, ledger

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT
    / "certificates"
    / "polynomial-vacuum-complete-two-loop-canonical-amplitude.json"
)
PARENT_SHA = "bdca7ba6d446a2afdc2abe197147e11be5e5d1d8c9bb118bd90867be9f947913"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_two_loop_canonical/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen complete two-loop light-pole report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_126_fully_rebuilt": PARENT_SHA,
        "same_actual_canonical_model_and_fixed_light_OS_conditions": True,
        "same_all_192_raw_graph_and_once_fixed_finite_contact_bounds": True,
        "same_entire_regulated_and_recursive_heavy_parameter_references": True,
        "canonical_not_derivative_coordinate_or_original_P8_completion": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = {k: bool(v) for k, v in audit.gates().items()}
    if not all(v is True for v in gates.values()):
        raise ValueError("A complete canonical two-loop amplitude gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.127.POLYNOMIAL_VACUUM_COMPLETE_TWO_LOOP_CANONICAL_AMPLITUDE",
        "date": "2026-09-10",
        "status": "COMPLETE_CANONICAL_TWO_LOOP_ELASTIC_FORWARD_COEFFICIENT_WITH_ALL_MARKED_CORE_COUNTERTERM_PRODUCTS_FIXED_LOCAL_REFERENCES_AND_LIGHT_LSZ; NOT_TWO_LOOP_SOURCE_MATCHING_ALL_HIGHER_LOOP_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/contractions.md",
            "notes/products.md",
            "notes/references.md",
            "notes/fields.md",
            "notes/assembly.md",
            "notes/bounds.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": len(rows),
        "proof_checks": gates,
        "complete_marked_graph_and_insertion_census": serialize(
            payload(contractions.data())
        ),
        "full_local_reference_and_cubic_product_ledger": serialize(
            payload(ledger.data())
        ),
        "canonical_LSZ_and_complete_order_assembly": serialize(
            {
                "fields": payload(fields.data()),
                "assembly": payload(assembly.data()),
            }
        ),
        "actual_complete_canonical_b2_bounds": serialize(
            {
                "actual": payload(calibration.data()),
                "tree": calibration.point(0),
                "one_loop": calibration.point(1),
                "two_loops": calibration.point(2),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete canonical two-loop elastic forward coefficient is bounded after an independent full external-label counterterm and LSZ ledger. All 1152 marked logarithmic vertex/heavy-mass cores in 4608 labelled graph records match 144 independently generated one-loop counterterm cographs, and all 72 disjoint-core pairs match 12 tree counterterm-product cographs, including signs. Fundamental cubic G bookkeeping retains the cubic-counterterm square exactly once. Inherited local references, the fixed finite quartic insertion, new constant potential contact and the complete light on-shell normalization account for every order-two contribution. The full canonical two-loop correction has absolute upper bound below 3 times 10^-607 and 10^-7 of tree. The combined one-plus-two-loop relative bound is below 10^-6, so the canonical forward coefficient through two loops is strictly positive. This is not an all-orders dispersion or UV-admissibility verdict, a two-loop source-aware derivative-map result or original P8 closure.",
        "not_established": [
            "The exact signed two-loop coefficient or all-higher-loop truncation error",
            "All-energy absorptive integrals complex contours or nonperturbative UV completion",
            "Two-loop source-aware derivative-coordinate matching or cosmological loop control",
            "Finite-gravity IR/Regge Delta absolute coupled renormalization or a common controlled bounce parent",
            "Full V/G/B classification any excluded original P8 row or original P8 closure",
        ],
        "verification_boundary": "Native exhaustive external-label marked-graph contraction and independent counterterm generation with exact signed Wick weights, actual local reference coefficients, complete cubic/mass product expansion, canonical bare-field and LSZ identities, disjoint raw-group ownership and actual rational error bounds. The Wick-contraction correspondence, local renormalization and functional even-vacuum/LSZ assembly are source-pinned written proofs, not proof-assistant formalized, nonperturbative or peer reviewed. Every ancestor, own source and report field is rebuilt read-only with unmodified scientific SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete canonical two-loop report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.127.POLYNOMIAL_VACUUM_COMPLETE_TWO_LOOP_CANONICAL_AMPLITUDE replay passed; full V/G/B and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
