"""Read-only complete all-finite tree coefficient bound certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_all_multiplicity_tree_source import verify as tree_input
from p8_vacuum_affine_uniform_soft_current_bound import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-uniform-all-tree-bound.json"
PARENT_SHA = "24258b5a7ce3e849edfcdea448a1f45a0b18bcf765f67c14b4cc7a45b19e84af"
TREE_SHA = "c0786edc77e68ff1c1a34b788674a782282ab96fb74be8268a345d5109840c40"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_uniform_all_tree_bound/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen uniform pure-soft current parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(tree_input.REPORT) != TREE_SHA:
        raise ValueError("The frozen all-multiplicity tree-source input changed")
    tree_input.validate_report(
        json.loads(tree_input.REPORT.read_text()), tree_input.build_report()
    )
    return {
        "S6_318_uniform_all_finite_soft_current_bound_rebuilt": PARENT_SHA,
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
        raise ValueError("A uniform complete-tree bound proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.319.UNIFORM_ALL_FINITE_MULTIPLICITY_COMPLETE_TREE_BOUND",
        "date": "2026-09-16",
        "status": "SCOPED_ALL_FINITE_COMPLETE_BARE_TREE_BOUND; NOT_INCLUSIVE_RATE_OR_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/core.md",
            "notes/vertices.md",
            "notes/majorant.md",
            "notes/calibration.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_hard_core": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_vertex_energy_and_complete_tree_majorant": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete original four-scalar/N-graviton finite tree is bounded relative to its positive unexpanded Born amplitude by3*n^2*N!*(2*10^16)^N/[kappa^(N/2)*product(w_i)]. The unique three-core decomposition, all-valence vertex bounds, paired hard propagators and disjoint soft-energy charges compose with the S318 current majorant. The estimate is uniform over angular and energy hierarchies on the generic physical domain, retaining explicit soft poles. Independent core inventories and the complete original5116-tree amplitude supplement the all-N proof. This is not an infrared-finite inclusive probability, quantum state, Regge estimate or original V/G/B/P8 closure.",
        "not_established": [
            "An infrared-finite inclusive all-N probability or convergent physical series",
            "The full overlapping all-N soft subtraction and real-virtual completion",
            "Complete finite hard real-virtual and evanescent matching",
            "An interacting quantum state, unitarity or absolute complex Regge",
            "The original common-parent bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Written all-N maximal-soft/core bijection; literal independent labeled core cuts through N5 and composition recovering full counts; all-valence vertex and hard-cut inequalities, exact trace-reversal norm identity and disjoint energy charges; nonnegative sequence/set majorants checked by independent rational coefficients; exact positive composition barrier and finite coefficient checks; complete original5116-tree bound and temporal gauge equivalence. Finite checks supplement the structural proof and do not establish infrared cancellation or a probability sum. Original SymPy is retained outside the FULL-only adapter.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The uniform complete-tree bound report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.319 uniform complete-tree bound replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
