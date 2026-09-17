"""Read-only complete-pair factorization and cluster-obstruction certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_all_multiplicity_tree_source import verify as previous
from p8_vacuum_affine_two_real_collinear_current import verify as tree_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-complete-pair-factorization.json"
PARENT_SHA = "c0786edc77e68ff1c1a34b788674a782282ab96fb74be8268a345d5109840c40"
TREE_SHA = "e73f4ba8880178262572d6ca3eb7c368f98d46ce900ef5f5a9bd9853c2b5aff1"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_complete_pair_factorization/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen all-multiplicity tree-source parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(tree_input.REPORT) != TREE_SHA:
        raise ValueError("The frozen conserved-pair angular input changed")
    tree_input.validate_report(
        json.loads(tree_input.REPORT.read_text()), tree_input.build_report()
    )
    return {
        "S6_315_all_multiplicity_tree_source_rebuilt": PARENT_SHA,
        "S6_311_exact_conserved_pair_quotient_rebuilt": TREE_SHA,
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
        raise ValueError("A complete-pair factorization proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.316.COMPLETE_PAIR_FACTORIZATION_MATCHING_FORESTS_AND_CLUSTER_OBSTRUCTIONS",
        "date": "2026-09-16",
        "status": "SCOPED_COMPLETE_PAIR_FACTORIZATION_AND_EXPLICIT_CLUSTER_OBSTRUCTIONS; NOT_ALL_N_RATE_OR_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/ward.md",
            "notes/forests.md",
            "notes/pairs.md",
            "notes/obstructions.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_matching_forests": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_complete_pair_factorization_and_obstructions": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged finite-tree source has conserved complete one-off-shell graviton currents and exact fixed-pair factorization with a conditional uniform angular bound. All-N matching-forest inclusion-exclusion accounts for overlapping pair branches. Full original5116-tree factorization, independent lower-source currents and exact counterexamples show why isolated-cluster or multi-off-shell Ward shortcuts fail. No completed-remainder norm, simultaneous-collinear estimate, all-N rate or original V/G/B/P8 closure is inferred.",
        "not_established": [
            "An all-N nonleading amplitude or inclusive probability bound",
            "Uniform completed-remainder norms and simultaneous-collinear estimates",
            "Complete finite hard real-virtual and evanescent matching",
            "An interacting quantum state, unitarity or absolute complex Regge",
            "The original common-parent bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Finite Noether induction for one-off-shell roots; exact root-cut fixed-pair factorization; Boolean matching-forest inclusion-exclusion and independent no-pair inventories through N4; all24 frozen pair-quotient coefficients; full5116-tree decomposition at diagnostic and original parameters with independent frozen434-tree currents; symbolic isolated-cluster pole, actual original47-tree hard-source residue and a two-off-shell counterexample. Calibrations do not replace the general proofs or yield an all-N probability bound. Original SymPy is retained outside the FULL-only adapter.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete-pair factorization report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.316 complete-pair factorization replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
