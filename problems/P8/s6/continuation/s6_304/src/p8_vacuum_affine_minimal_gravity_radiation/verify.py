"""Read-only full Einstein radiation and uniform fixed-resolution real-rate certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_finite_forward_phase import verify as previous
from p8_vacuum_affine_one_newton_inclusive_assembly import verify as born_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-minimal-gravity-radiation.json"
PARENT_SHA = "d02d6507a596c20adad3d3fd5a654f1a8abdca1857a1a9b49bc4c91bc5f2457a"
BORN_SHA = "6206f5d56959f01a8467e2c2b74ff43a2fad3cc75df0cee67e9bb8ce82f39564"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_minimal_gravity_radiation/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen finite forward-phase parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(born_input.REPORT) != BORN_SHA:
        raise ValueError("The frozen full positive Born normalization changed")
    born_input.validate_report(
        json.loads(born_input.REPORT.read_text()), born_input.build_report()
    )
    return {
        "S6_303_whole_finite_forward_phase_and_known_rate_rebuilt": PARENT_SHA,
        "S6_297_full_positive_Born_normalization_rebuilt": BORN_SHA,
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
        raise ValueError("A gravity-radiation proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.304.COMPLETE_EINSTEIN_RADIATION_GENERAL_WARD_IDENTITY_AND_UNIFORM_FIXED_RESOLUTION_FULL_TREE_ERROR",
        "date": "2026-09-15",
        "status": "SCOPED_FULL_SELECTED_REAL_TREE_AND_UNIFORM_FIXED_RESOLUTION_ERROR; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/vertices.md",
            "notes/ward.md",
            "notes/bounds.md",
            "notes/matching.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_canonical_vertices": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_general_Ward_identity_and_real_rate": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "All21 minimal Einstein radiation graphs obey a general on-shell channel Ward identity and match the canonical Born soft limit. Together with the26 matter graphs, the full47-graph selected real tree has a finite real-minus-soft normalized rate error below1e32/kappa=1e-768 for every fixed0<x<=1/8, at every nonforward hard angle. Finite forward cross section, finite hard matching, quantum all-N/all-loop/Regge and original V/G/B/P8 remain OPEN.",
        "not_established": [
            "Values or bounds of the three physical finite matching coordinates",
            "Full detector-inclusive hard rate beyond the computed selected corrections",
            "Finite forward cross section or loop-square and all-order control",
            "Complex transfer analyticity, high-energy Regge and original quantum construction",
            "Original V/G/B/P8 closure",
        ],
        "verification_boundary": "Canonical action expansion, general invariant-index Ward polynomial, independent exact component and metric-expansion calibrations, and written propagator/component/phase-space estimates. The full selected tree and its interferences are retained. A physical recoil transfer gap and a two-region energy estimate control the normalized real correction at fixed resolution and all nonforward hard angles, without asserting a finite forward cross section. Native/direct/ordinary/CLI use original SymPy, the adapter is FULL-only, frozen ancestors and original statuses are unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The gravity-radiation report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.304 full Einstein radiation replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
