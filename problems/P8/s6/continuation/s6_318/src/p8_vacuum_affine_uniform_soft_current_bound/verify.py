"""Read-only all-finite-multiplicity weighted temporal soft-current certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_all_multiplicity_tree_source import verify as tree_input
from p8_vacuum_affine_temporal_tree_reorganization import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-uniform-soft-current-bound.json"
PARENT_SHA = "477509f3e37369b43fc2abdf0db75a498fa222b476ad23ed56ca0f6fbbe3b5fd"
TREE_SHA = "c0786edc77e68ff1c1a34b788674a782282ab96fb74be8268a345d5109840c40"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_uniform_soft_current_bound/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen gauge-complete temporal-tree parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(tree_input.REPORT) != TREE_SHA:
        raise ValueError("The frozen all-multiplicity tree-source input changed")
    tree_input.validate_report(
        json.loads(tree_input.REPORT.read_text()), tree_input.build_report()
    )
    return {
        "S6_317_complete_temporal_tree_reorganization_rebuilt": PARENT_SHA,
        "S6_315_unchanged_all_multiplicity_tree_source_rebuilt": TREE_SHA,
        "S279_remains_rejected_and_archived": True,
        "S275_S276_S277_and_scoped_P8a_unchanged": True,
        "all6_historical_physical_qualifications_unchanged": True,
        "all_original_primitive_and_matching_frontiers_retained": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A uniform soft-current bound proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.318.ALL_FINITE_MULTIPLICITY_WEIGHTED_TEMPORAL_SOFT_CURRENT_BOUND",
        "date": "2026-09-16",
        "status": "SCOPED_ALL_FINITE_WEIGHTED_TEMPORAL_SOFT_CURRENT_BOUND; NOT_INCLUSIVE_RATE_OR_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/geometry.md",
            "notes/grading.md",
            "notes/majorant.md",
            "notes/calibration.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_weighted_geometry": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_null_grading_and_all_finite_current_majorant": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged complete temporal pure-soft current has a weighted norm bound uniform over all angular and nested collinear configurations at every finite multiplicity, with explicit W^n/product(w_i) energy dependence. Arbitrary Rosen geometry and reflection grading remove the dangerous root powers; the conserved temporal inverse closes the weighted norm. A labeled recurrence and exact counting-series barrier give C_n<=2*(2*10^10)^(n-1)*n!. Complete26-tree collinear families and low-valence failed controls test the implementation. No complete hard amplitude, all-N inclusive probability, quantum state, Regge or original V/G/B/P8 closure is inferred.",
        "not_established": [
            "A complete hard four-scalar/all-N nonleading amplitude or inclusive probability bound",
            "The full overlapping all-N soft subtraction and real-virtual completion",
            "Complete finite hard real-virtual and evanescent matching",
            "An interacting quantum state, unitarity or absolute complex Regge",
            "The original common-parent bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Independent arbitrary Rosen Ricci/Gamma-Gamma calculation; exact tensor and momentum grading; EH3/EH4/EH5 polynomial calibrations and a nonzero wrong-weight control; written coefficient-l1 Banach transfer of the frozen all-order vertex bound; complete conserved-root weighted inverse; labeled-partition versus EGF counts through N10 and an exact positive barrier; two complete26-tree collinear families testing the weighted field, velocity and longitudinal components at both endpoints. Finite checks do not replace the all-order structural proof or establish an inclusive rate. Original SymPy is retained outside the FULL-only adapter.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The uniform soft-current bound report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.318 uniform soft-current bound replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
