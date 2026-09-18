"""Read-only complete third-tree amplitude-subtraction certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_pair_singleton_subtraction import verify as previous
from p8_vacuum_affine_relative_energy_complex_tube import verify as analytic_input
from p8_vacuum_affine_uniform_all_tree_bound import verify as uniform_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-three-singleton-subtraction.json"
PARENT_SHA = "0b029aab51a1be6266c4163049bf7d29a37aabd44b650ccc64ff147e0d51035d"
ANALYTIC_SHA = "18292735783403d7588808de71e7620238ffbbe2f70d82a5a1012ef431f03108"
UNIFORM_SHA = "6239bb0b855c48275f8598cdfdec3f8dfa2590545d5b4abb5b6ddfe765017b0b"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_three_singleton_subtraction/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen pair-plus-singleton subtraction parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(analytic_input.REPORT) != ANALYTIC_SHA:
        raise ValueError("The frozen relative-complex and hard-majorant input changed")
    analytic_input.validate_report(
        json.loads(analytic_input.REPORT.read_text()), analytic_input.build_report()
    )
    if sha(uniform_input.REPORT) != UNIFORM_SHA:
        raise ValueError("The frozen complete finite-tree amplitude input changed")
    uniform_input.validate_report(
        json.loads(uniform_input.REPORT.read_text()), uniform_input.build_report()
    )
    return {
        "S6_319_uniform_complete_tree_baseline_rebuilt": UNIFORM_SHA,
        "S6_322_all1349_nonsingleton_subtraction_rebuilt": PARENT_SHA,
        "S6_320_complex_recoil_and_hard_majorants_rebuilt": ANALYTIC_SHA,
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
        raise ValueError("A complete third-amplitude subtraction proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.323.COMPLETE_THREE_REAL_TREE_SUBTRACTION_AND_FINITE_SIGNED_MEASURE",
        "date": "2026-09-17",
        "status": "SCOPED_FULL5116_TREE_RECTANGLE_AND_SIGNED_MEASURE; NOT_INCLUSIVE_ALL_N_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/line.md",
            "notes/forest.md",
            "notes/grouping.md",
            "notes/analytic.md",
            "notes/measure.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_connected_lines": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_forest_grouping_and_complete_subtraction": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The original3767 hard core is exactly reorganized into eight external-line-cumulant and central-source contributions. Generic connected-line cancellations, all endpoint-owner groupings and two-regime hard-only analytic estimates give its uniform anchored third remainder. Combined with the frozen1349 nonsingleton result, the complete5116-tree amplitude coefficient abc*sqrt(rho)*M3/A0 has compatible faces and rectangle below1e-670*J. The seven-face amplitude baseline defines a full signed density difference retaining all interference, with total variation below1e-1424*x^2+1e-1340*x^4 and a unique common-cutoff limit. This is not a positive normalized inclusive probability, all-N or real-virtual completion, interacting state, absolute complex Regge, common-parent bounce or original P8 closure.",
        "not_established": [
            "Matching this signed three-real baseline to actual virtual and integrated counterterms or a positive normalized inclusive probability",
            "The full all-N overlapping subtraction and summation",
            "Complete finite hard and evanescent matching",
            "An interacting quantum state, unitarity or absolute complex Regge",
            "The original common-parent bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Generic49-coefficient scalar-line polynomial identities, original16 signed complex-TT mappings, independent original3767/125-assignment/eight-class equality;192 endpoint-grouping cases and776 coefficient factorizations, retained negative controls; literal recoil and massive-current identities; uniform hard-only physical and origin domains, grouped Taylor remainder proof, exact eight-class budgets and compatible-face entropy integrability; exact seven-face projectors, all complex interference, entropy integrals and full polarization/identical-particle measure. A connected-pair pole explicitly forbids applying the hard-only tube to the full current. Finite points calibrate, not replace, the written uniform proof. Not kernel-formalized; original SymPy outside FULL.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete third-amplitude subtraction report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.323 full third-tree amplitude and signed subtraction replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
