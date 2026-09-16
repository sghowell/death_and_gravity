"""Read-only off-shell current and uniform selected two-real collinear bound."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_complete_two_graviton_tree import verify as previous
from p8_vacuum_affine_minimal_gravity_radiation import verify as radiative_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-two-real-collinear-current.json"
PARENT_SHA = "e60342ad8eacf8ef94a07d510fe12cfcbf71926fb0977d0aeb0b794b31154732"
RADIATIVE_SHA = "a34d5160c55891304eda5f83bf90970854774e8aae327d9375b89e5fa797063e"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_two_real_collinear_current/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen complete two-graviton tree parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(radiative_input.REPORT) != RADIATIVE_SHA:
        raise ValueError("The frozen canonical Einstein radiation input changed")
    radiative_input.validate_report(
        json.loads(radiative_input.REPORT.read_text()), radiative_input.build_report()
    )
    return {
        "S6_310_complete_two_graviton_tree_rebuilt": PARENT_SHA,
        "S6_304_canonical_Einstein_radiation_rebuilt": RADIATIVE_SHA,
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
        raise ValueError("An offshell-current/collinear proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.311.COMPLETE_OFFSHELL_CURRENT_AND_UNIFORM_TWO_REAL_COLLINEAR_SECTOR_BOUND",
        "date": "2026-09-16",
        "status": "SCOPED_UNIFORM_SELECTED_TWO_REAL_COLLINEAR_CURRENT; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/ward.md",
            "notes/factorization.md",
            "notes/collinear.md",
            "notes/bounds.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_offshell_current": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_collinear_factorization_and_uniform_bound": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged full434-tree amplitude factors into387 complementary graphs and47 graphs carrying the q1+q2 propagator. The latter use a generally conserved complete timelike-root current. Exact physical polarization contraction cancels its apparent collinear angular pole and gives a uniform coefficient budget. A conservative original full-Born current estimate yields|Mpair|/(Am+AG)<9e-400(1/a+1/b) uniformly in emitted and nonforward hard directions. Soft overlaps and the complementary amplitude are not omitted-error bounded here; integrated two-real/all-N rates, matching, Regge and original V/G/B/P8 remain OPEN.",
        "not_established": [
            "A bound on the complementary387 graphs or full434-tree remainder",
            "An IR-finite integrated two-real correction, overlap subtraction or virtual pairing",
            "All-N detector rates, radiative hard loops and finite physical matching",
            "Unitarity, complex Regge control or the constructed original quantum state",
            "The common-parent bounce or original V/G/B/P8 closure",
        ],
        "verification_boundary": "General off-shell invariant-index Ward polynomial, independent26-matter tensor, exact pair-current factorization at original parameters, all24 physical pair coefficients and their polynomial absolute-sum bounds, explicit recoil/propagator/component budgets and exact extreme-angle calibrations. These prove a selected-class amplitude bound, not integrated infrared finiteness or all-order physics. Native/direct/ordinary/CLI retain original SymPy; the adapter is FULL-only. All frozen ancestors and scoped P8(a) remain unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The two-real collinear current report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.311 two-real collinear current replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
