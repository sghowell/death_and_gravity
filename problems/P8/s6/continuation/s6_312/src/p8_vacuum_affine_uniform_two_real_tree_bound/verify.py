"""Read-only original complete two-real tree bound; not IR-finite P8 closure."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_complete_two_graviton_tree import verify as tree_input
from p8_vacuum_affine_two_real_collinear_current import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-uniform-two-real-tree-bound.json"
PARENT_SHA = "e73f4ba8880178262572d6ca3eb7c368f98d46ce900ef5f5a9bd9853c2b5aff1"
TREE_SHA = "e60342ad8eacf8ef94a07d510fe12cfcbf71926fb0977d0aeb0b794b31154732"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_uniform_two_real_tree_bound/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen two-real collinear-current parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(tree_input.REPORT) != TREE_SHA:
        raise ValueError(
            "The frozen canonical complete two-graviton tree input changed"
        )
    tree_input.validate_report(
        json.loads(tree_input.REPORT.read_text()), tree_input.build_report()
    )
    return {
        "S6_311_two_real_collinear_current_rebuilt": PARENT_SHA,
        "S6_310_complete_two_graviton_tree_rebuilt": TREE_SHA,
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
        raise ValueError("A complete two-real tree bound proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.312.UNIFORM_ORIGINAL_COMPLETE_TWO_REAL_TREE_AMPLITUDE_BOUND",
        "date": "2026-09-16",
        "status": "SCOPED_UNIFORM_COMPLETE434_BARE_TREE_BOUND; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/topology.md",
            "notes/propagators.md",
            "notes/vertices.md",
            "notes/bounds.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_regular_topology": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_common_gaps_and_complete_tree_bound": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The original434-tree amplitude has a uniform upper bound |M6|/(Am+AG)<2e-395/(ab) for E in[5/4,2],a,b>0,a+b<=1/8, all physical emitted directions and all nonforward hard Born angles. Every one of the387 regular graphs is covered by common physical cluster gaps and canonical scalar/Einstein vertex budgets; the frozen47-graph collinear-current theorem controls the remaining pair class. This is a bare soft-divergent amplitude envelope, not an infrared-finite integrated rate. Soft overlaps, virtual pairing, all-N rates, finite matching, Regge and original V/G/B/P8 remain OPEN.",
        "not_established": [
            "A separately physical detector observable for either gauge-fixed graph class",
            "An IR-finite integrated two-real correction, overlap subtraction or virtual pairing",
            "All-N detector rates, radiative hard loops and finite physical matching",
            "Unitarity, complex Regge control or the constructed original quantum state",
            "The common-parent bounce or original V/G/B/P8 closure",
        ],
        "verification_boundary": "General Kallen/energy-gradient/recoil inequalities for arbitrary radiation assignments; exact inventory of every434 labeled tree and387 regular denominator profiles; canonical multilinear vertex coefficient budgets; exact original-parameter complete434 sums at extreme angular boundaries. Samples calibrate the written uniform proof rather than replace it. Squared soft-envelope divergence is retained explicitly. Native/direct/ordinary/CLI retain original SymPy; the adapter is FULL-only. Frozen ancestors and scoped P8(a) are unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete two-real tree bound report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.312 uniform complete two-real tree bound replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
