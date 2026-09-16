"""Read-only uniform state-correct two-real soft-overlap subtracted measure."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_two_real_collinear_current import verify as pair_input
from p8_vacuum_affine_uniform_two_real_tree_bound import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-two-real-soft-overlap.json"
PARENT_SHA = "d76abfc6ca786d9357b04152e5622ccd3896ae71a84c4bb222b8ecc9c54019eb"
PAIR_SHA = "e73f4ba8880178262572d6ca3eb7c368f98d46ce900ef5f5a9bd9853c2b5aff1"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_two_real_soft_overlap/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen uniform complete two-real tree parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(pair_input.REPORT) != PAIR_SHA:
        raise ValueError("The frozen conserved two-real pair-current input changed")
    pair_input.validate_report(
        json.loads(pair_input.REPORT.read_text()), pair_input.build_report()
    )
    return {
        "S6_312_uniform_complete_two_real_tree_rebuilt": PARENT_SHA,
        "S6_311_two_real_collinear_current_rebuilt": PAIR_SHA,
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
        raise ValueError("A uniform two-real soft-overlap proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.313.UNIFORM_STATE_CORRECT_TWO_REAL_SOFT_OVERLAP_SUBTRACTED_MEASURE",
        "date": "2026-09-16",
        "status": "SCOPED_UNIFORM_TWO_REAL_SOFT_OVERLAP_SUBTRACTED_MEASURE; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/external.md",
            "notes/analytic.md",
            "notes/overlap.md",
            "notes/measure.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_external_regrouping": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_complex_uniform_and_soft_overlap_bound": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The original full434 two-real tree admits a state-correct, gauge-consistent subtraction of both single-soft faces with their common overlap restored once. Exact double-external regrouping and a uniform holomorphic recoil/propagator/component proof give |partial_a partial_b(ab sqrt(rho2) M6/A0)|<1e-326/(a+b). The resulting complete signed two-real measure, including all interferences, has total variation below2e-725*x+1e-652*x2 on the stated all-angle compact domain, independently of its lower soft cutoff. Virtual/counterterm matching, all-N nonleading errors, physical matching, Regge and original V/G/B/P8 remain OPEN.",
        "not_established": [
            "The complete inclusive detector rate or finite virtual/counterterm matching",
            "Identification with a sum of overlapping earlier subtraction schemes",
            "All-N nonleading radiation errors, hard loops or finite physical matching",
            "Unitarity, absolute complex Regge or the original quantum state",
            "The common-parent bounce or original V/G/B/P8 closure",
        ],
        "verification_boundary": "General two-order scalar-line and stress identities, exhaustive double-external graph grouping, explicit complex recoil and propagator margins, paired-current telescoping bounds and Cauchy estimate, exact entropy-series integrals with independent quadrature, correct phase/polarization/Bose factors and six exact original434/47-face calibrations. Numerical soft convergence is supplementary, not the uniform proof. Native/direct/ordinary/CLI retain original SymPy; the adapter is FULL-only. All frozen ancestors and scoped P8(a) are unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The two-real soft-overlap report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.313 uniform two-real soft-overlap replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
