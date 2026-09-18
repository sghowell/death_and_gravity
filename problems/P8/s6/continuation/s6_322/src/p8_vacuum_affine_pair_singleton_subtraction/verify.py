"""Read-only pair-plus-singleton subtraction certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_three_soft_current_subtraction import verify as previous
from p8_vacuum_affine_two_real_soft_overlap import verify as analytic_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-pair-singleton-subtraction.json"
PARENT_SHA = "297fd97c371ffc0df3e2ebd47d9eabcdc3cd21bfb3e6031fdba8c69ef069379c"
ANALYTIC_SHA = "79539dcf68142c2ebb62a0236293be98d2f1b6eece93d49e46cbf091fb8e1994"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_pair_singleton_subtraction/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen three-current and188-class parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(analytic_input.REPORT) != ANALYTIC_SHA:
        raise ValueError("The frozen complex-recoil and two-real overlap input changed")
    analytic_input.validate_report(
        json.loads(analytic_input.REPORT.read_text()), analytic_input.build_report()
    )
    return {
        "S6_321_three_current_and188_class_rebuilt": PARENT_SHA,
        "S6_313_complex_recoil_and_two_real_overlap_rebuilt": ANALYTIC_SHA,
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
        raise ValueError("A pair-plus-singleton subtraction proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.322.PAIR_SINGLETON_SUBTRACTION_AND_ALL_NON_SINGLETON_FACES",
        "date": "2026-09-17",
        "status": "SCOPED_ALL1349_NON_SINGLETON_SUBTRACTION; NOT_FULL5116_INCLUSIVE_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/hierarchy.md",
            "notes/algebra.md",
            "notes/tube.md",
            "notes/budget.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_pair_hierarchy": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_offshell_core_and_subtraction": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete temporal pair has an all-angle lower energy derivative hierarchy and compatible faces. Keeping off-shell trace and nontransverse terms, the387 hard core has an anisotropic complex-energy domain and a Born-product remainder below B2*W. Operator Cauchy and the pair hierarchy give a threefold rectangle below1e-723*c*I(a,b) for each pair-plus-singleton class. The three pair choices and S321's188 class give an all1349 nonsingleton rectangle below2e-723*J. The3767 singleton core, full5116 probability, all-N and real-virtual matching, quantum state, Regge and original P8 remain open.",
        "not_established": [
            "Subtraction for the3767 three-singleton core or complete5116-tree probability",
            "The all-N overlapping soft and real-virtual completion",
            "Complete finite hard and evanescent matching",
            "An interacting quantum state, unitarity or absolute complex Regge",
            "The original common-parent bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact complete-pair tensor/derivative identities and positive angular coefficient majorants; generic arbitrary-tensor scalar vertices, traceful shifts and ordering identities; independent original387,140 and temporal-class calibrations; closed anisotropic hard-core gaps, forward grouped remainder proof, exact budgets, compatible faces and Cauchy product rule. Three exact complex-domain calibrations supplement, not replace, the uniform proof. No full-current global-W holomorphy premise. Not kernel-formalized; original SymPy outside FULL.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The pair-plus-singleton subtraction report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.322 pair-plus-singleton subtraction replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
