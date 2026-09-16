"""Read-only complete selected two-graviton formal tree certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_minimal_gravity_radiation import verify as radiative_input
from p8_vacuum_affine_single_residual_soft_dressing import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-complete-two-graviton-tree.json"
PARENT_SHA = "22e89710b44e3df51ab35cf4ef541f70c42c68ee740a29f261f9f9d16694b59c"
RADIATIVE_SHA = "a34d5160c55891304eda5f83bf90970854774e8aae327d9375b89e5fa797063e"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_complete_two_graviton_tree/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen single-residual dressing parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(radiative_input.REPORT) != RADIATIVE_SHA:
        raise ValueError("The frozen complete one-graviton radiation input changed")
    radiative_input.validate_report(
        json.loads(radiative_input.REPORT.read_text()), radiative_input.build_report()
    )
    return {
        "S6_309_complete_state_correct_single_residual_dressing_rebuilt": PARENT_SHA,
        "S6_304_complete_Einstein_radiation_rebuilt": RADIATIVE_SHA,
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
        raise ValueError("A complete two-graviton tree proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.310.COMPLETE_SELECTED_TWO_GRAVITON_TREE_AND_STATE_CORRECT_LEADING_SOFT_LIMITS",
        "date": "2026-09-16",
        "status": "SCOPED_COMPLETE_SELECTED_TWO_REAL_FORMAL_TREE; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/vertices.md",
            "notes/recursion.md",
            "notes/ward.md",
            "notes/soft.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_rooted_trees": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_action_and_soft_calibrations": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged selected four-scalar/two-graviton formal tree is completely constructed from the covariant metric jets by a finite rooted-tree recursion, with434 independently counted graphs. Literal higher-vertex action checks and frozen lower-point calibration fix normalization. Complete six-point Ward/Bose and original-parameter simultaneous/hierarchical soft checks pass. Fixed-generic-kinematics tree soft limits use the full marked radiative state. No uniform integrated two-real remainder, all-N detector rate, physical matching, Regge or original V/G/B/P8 closure is asserted.",
        "not_established": [
            "A uniform integrated two-real nonleading residual and its overlapping soft subtraction",
            "Uniform collinear/forward limits or virtual pairing for the complete two-real correction",
            "The full all-N physical detector rate, hard radiative loops and finite matching",
            "Unitarity, high-energy complex Regge control or a constructed quantum state",
            "The common-parent bounce or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact canonical metric-jet derivation, root-cut bijection, independent formal graph enumeration, full frozen lower-point comparisons, non-opposite six-point Ward/Bose tests, wrong-vertex controls and exact original-parameter soft calibrations. General tree Ward decoupling and fixed-kinematics soft statements use the written covariant-action arguments; finite samples do not prove uniform error estimates. Native/direct/ordinary/CLI keep original SymPy; the adapter is FULL-only. Frozen ancestors and scoped P8(a) remain unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete two-graviton tree report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.310 complete two-graviton tree replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
