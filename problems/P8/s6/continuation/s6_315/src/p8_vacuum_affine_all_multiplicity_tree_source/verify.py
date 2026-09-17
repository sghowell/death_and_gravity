"""Read-only arbitrary finite-multiplicity canonical tree-source certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_complete_two_graviton_tree import verify as tree_input
from p8_vacuum_affine_state_transport_soft_matching import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-all-multiplicity-tree-source.json"
)
PARENT_SHA = "ea01ed8699689c68db6d52ed9c9ebf25ee0015e478480e43391a60a68a03dd56"
TREE_SHA = "e60342ad8eacf8ef94a07d510fe12cfcbf71926fb0977d0aeb0b794b31154732"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_all_multiplicity_tree_source/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen state-transport matching parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(tree_input.REPORT) != TREE_SHA:
        raise ValueError("The frozen lower-multiplicity canonical tree input changed")
    tree_input.validate_report(
        json.loads(tree_input.REPORT.read_text()), tree_input.build_report()
    )
    return {
        "S6_314_state_transport_matched_soft_reference_rebuilt": PARENT_SHA,
        "S6_310_complete_lower_multiplicity_tree_source_rebuilt": TREE_SHA,
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
        raise ValueError("An arbitrary-multiplicity tree-source proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.315.ARBITRARY_FINITE_MULTIPLICITY_CANONICAL_TREE_SOURCE_AND_MAJORANTS",
        "date": "2026-09-16",
        "status": "SCOPED_ALL_FINITE_TREE_SOURCE_AND_VERTEX_TOPOLOGY_MAJORANTS; NOT_ALL_N_RATE_OR_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/vertices.md",
            "notes/topology.md",
            "notes/trees.md",
            "notes/bounds.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_all_order_vertices": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_all_N_topology_and_finite_tree_source": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged selected covariant action has an exact canonical vertex and labeled rooted-tree construction for every finite number of real gravitons. A closed all-N graph count, factorial-exponential vertex majorants and an explicit graph-count majorant are proved. Independent literal EH5/scalar-fourth action checks and a complete original5116-tree three-real amplitude supplement the proof. These separate majorants do not control the full all-N nonleading radiation sum, hard loops, physical matching or original V/G/B/P8.",
        "not_established": [
            "An all-N nonleading amplitude or inclusive probability bound",
            "Quantitative arbitrary soft-subtree collinear cancellation estimates",
            "Complete finite hard real-virtual and evanescent matching",
            "An interacting quantum state, unitarity or absolute complex Regge",
            "The original common-parent bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "General mixed-coefficient action derivation, finite root-cut bijection and Noether Ward induction; a solved formal graph generating function and analytic coefficient majorant; separate norm bounds for all canonical vertices; literal determinant/adjugate/indexed EH5 checks; independent root-cut counts through N4; exact agreement with the frozen lower amplitudes; two full N3 gauge directions, nonzero higher-vertex omission controls and an original-parameter N3 amplitude. Finite calibrations are not the general proof or an all-N probability estimate. Original SymPy is retained outside the FULL-only adapter.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The arbitrary-multiplicity tree-source report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.315 arbitrary-multiplicity tree-source replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
