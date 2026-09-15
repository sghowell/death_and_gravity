"""Read-only minimal-gravity pole completion and massive soft-kernel report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_ordinary_gravity_ward_bridge import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-massive-gravity-pole-completion.json"
)
PARENT_SHA = "9bbef0eecaf83290e7e6fdae4d3dc302c9e21177b3ad57e87195f927165f876e"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_massive_gravity_pole_completion/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen ordinary Ward-bridge report changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_287_whole_ordinary_Ward_charge_continuity_and_entire_ancestry_rebuilt": PARENT_SHA,
        "S279_remains_rejected_and_archived": True,
        "S275_S276_S277_and_scoped_P8a_unchanged": True,
        "all6_historical_physical_qualifications_unchanged": True,
        "all_original_primitive_and_matching_frontiers_retained": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "A pole-completion, whole-graph or soft-kernel proof gate failed"
        )
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.288.ORIGINAL_MINIMAL_GRAVITY_ONE_LOOP_PHYSICAL_POLE_COMPLETION_MODULO_ONE_NEWTON_AND_TWO_LOCAL_ANCHORS_WITH_WHOLE_MASSIVE_COULOMB_KERNEL",
        "date": "2026-09-15",
        "status": "SCOPED_MINIMAL_GR_POLE_COMPLETION_AND_CONDITIONAL_MASSIVE_SOFT_KERNEL; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/poles.md",
            "notes/graphs.md",
            "notes/soft.md",
            "notes/matching.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_physical_pole_jets": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_graph_completion_and_retained_Coulomb_kernel": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The full-D cut and complete evanescent-tree LSZ physical principal parts have an exact four-function crossing completion. The entire tensor-box numerator reduces to its fixed scalar box plus propagator-cancelled remainders, and all minimal graph classes satisfy the required transfer-weighted limits. Whole-component ordinary Ward/LSZ continuity, the fixed metric volume density and retained Newton coefficient determine the missing physical pole terms, modulo two regular local coefficients. No finite anchor is chosen. The exact massive soft kernel retains a nonzero forward Coulomb phase and a strictly signed interior attenuation density; the scalar boxes remain. This is not the full-source finite amplitude, detector/Regge control or original V/G/B/P8 closure.",
        "not_established": [
            "Finite values or signs of the light Newton and two on-shell regular local matching coefficients",
            "Full-source crossed amplitude, improved b20, other C/g/heavy/Proca/M1 sectors or higher loops",
            "Finite F2 or F1 slope, an exact interacting LSZ pole, finite physical detector/dressing/hard-unitarity errors",
            "A finite-G forward-limit theorem, fixed-transfer Regge contour bound or physical EFT cutoff",
            "Original quantum state/domain/measure/bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact literal all-D full-source/master/leg/tensor identities, independent D4/D5/D6 current components, complete valence inventory, crossing polynomial and tagged UV/IR algebra support the written diagram/continuity/cut-completeness arguments. Exact parameter-integrand signs and sheet identities support the massive soft result. These analytic proofs are not FORMALIZED. Native/direct/ordinary/CLI retain original SymPy; exact-GCD is FULL-only. No frozen source/raw or finite prescription is changed.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The minimal-gravity pole completion report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.288 minimal gravity pole completion and massive Coulomb replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
