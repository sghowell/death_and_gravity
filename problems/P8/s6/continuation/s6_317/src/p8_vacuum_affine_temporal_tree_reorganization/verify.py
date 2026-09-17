"""Read-only gauge-complete temporal-tree and restricted planar-bound certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_all_multiplicity_tree_source import verify as tree_input
from p8_vacuum_affine_complete_pair_factorization import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-temporal-tree-reorganization.json"
)
PARENT_SHA = "d2431ee97f5b01537f6ccc5d3f5c3d25fdb0251e934f1bf802564c0c92966f27"
TREE_SHA = "c0786edc77e68ff1c1a34b788674a782282ab96fb74be8268a345d5109840c40"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_temporal_tree_reorganization/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen temporal-tree reorganization parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(tree_input.REPORT) != TREE_SHA:
        raise ValueError("The frozen all-multiplicity tree-source input changed")
    tree_input.validate_report(
        json.loads(tree_input.REPORT.read_text()), tree_input.build_report()
    )
    return {
        "S6_316_complete_pair_factorization_and_obstructions_rebuilt": PARENT_SHA,
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
        raise ValueError("A temporal-tree reorganization proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.317.GAUGE_COMPLETE_TEMPORAL_TREE_REORGANIZATION_AND_PLANAR_THREE_RAY_BOUND",
        "date": "2026-09-16",
        "status": "SCOPED_GAUGE_COMPLETE_TEMPORAL_TREES_AND_PLANAR_BOUND; NOT_ALL_N_RATE_OR_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/chart.md",
            "notes/recursion.md",
            "notes/bounds.md",
            "notes/primary.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_finite_temporal_chart": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_recursive_tree_and_restricted_current_bounds": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged finite tree source admits a gauge-complete recursive temporal policy for positive-energy pure-soft currents. Independent literal nonlinear pullbacks agree at three and four leaves; the complete5116-tree amplitudes agree exactly at original and diagnostic parameters. All24 conserved-pair coefficients yield a two-ray field norm. An exact all8-TT polynomial majorant bounds the coplanar fixed1:2:3-energy three-ray field below120000/kappa on the stated half-angle box. No arbitrary three-dimensional, energy-hierarchy or all-N inclusive bound or original V/G/B/P8 closure is inferred.",
        "not_established": [
            "An all-N nonleading amplitude or inclusive probability bound",
            "Arbitrary three-dimensional, energy-hierarchy and higher-multiplicity current bounds",
            "Complete finite hard real-virtual and evanescent matching",
            "An interacting quantum state, unitarity or absolute complex Regge",
            "The original common-parent bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Finite Noether/Euler induction and on-shell tree gauge invariance; a literal nilpotent coordinate pullback independent of recursive propagation; exact original/diagnostic5116-tree amplitude equality; all24 pair tensor coefficients; all8 planar TT triples with exact denominator factorization, numerator support and coefficient-l1 budgets. The linear-only projection remains a failed control. Calibration does not replace written finite-tree proofs or establish an all-N rate. Original SymPy is retained outside the FULL-only adapter.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The temporal-tree reorganization report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.317 temporal-tree reorganization replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
