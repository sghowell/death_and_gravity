"""Read-only original massive physical two-graviton cut report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_massive_detector_ir_cut import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-massive-graviton-pair-cut.json"
PARENT_SHA = "2cb6b7d65b559dcfc45dbaed1a8f058088cc15af87da32a9254ab70d9ca9afba"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_massive_graviton_pair_cut/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen original M1 cut and detector report changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_278_full_original_M1_cut_detector_and_entire_ancestry_rebuilt": PARENT_SHA,
        "S275_corrected_finite_hybrid_not_refuted": True,
        "S276_S277_paths_not_reidentified_as_the_original_state": True,
        "S278_M1_cut_and_conditional_soft_dictionary_unchanged": True,
        "all6_historical_physical_qualifications_unchanged": True,
        "all_original_primitive_and_matching_frontiers_retained": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A massive physical graviton cut proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.279.ORIGINAL_MASSIVE_TWO_GRAVITON_CUT_FROM_LITERAL_EINSTEIN_AND_MATTER_VERTICES_WITH_COMPLETE_ANGULAR_SEWING_ALL_EVEN_SPIN_SQUARES_AND_PHYSICAL_WINDOW_BOUND",
        "date": "2026-09-15",
        "status": "SCOPED_FORMAL_PHYSICAL_TWO_GRAVITON_CUT; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/amplitude.md",
            "notes/angular.md",
            "notes/spins.md",
            "notes/window.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_literal_graviton_pair_tree": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_physical_angular_spectral_and_window_cut": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The entire original fixed source admits explicit formal-loop bookkeeping without deleting its three vacuum constants. Literal EH cubic and scalar contact/exchange calculations establish all four canonically normalized two-graviton helicity amplitudes. The complete phase-sensitive physical cut has an exact massive forward expression, all-even-spin positive square weights and a uniform finite physical-window bound. The latter is neither b20 nor a complete low-cut subtraction or Regge allowance. The unphysical massless-threshold continuation, other quantum sectors, finite matching and original state/domain problem remain unresolved. Original V/G/B/P8 remain OPEN.",
        "not_established": [
            "An exact interacting quantum vacuum, physical LSZ mass or resummed finite-gravity amplitude",
            "The unphysical0<s<4m^2 cut, complete crossed low-cut subtraction or full b20",
            "All other species, real local matching, higher loops or finite detector and Regge errors",
            "A physical cutoff inferred from a formal integration cap, or original quantum-bounce/state completion",
            "Original V/G/B/P8 closure, all-parent exclusion or formal verification of written harmonic analysis",
        ],
        "verification_boundary": "Exact whole-source loop bookkeeping, literal Einstein and scalar vertices, both Ward identities, physical TT projection, full angular phase sewing and exact massive forward integration support the written all-spin and window arguments. Finite Q checks are not an exhaustive spin proof. The normalization is independently fixed from the source, not a printed spinor prefactor. Native/direct/ordinary/CLI retain original SymPy; only a captured complete FULL replay uses the audited exact-GCD adapter. All frozen science and qualifications remain unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete massive graviton pair cut report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.279 original massive physical two-graviton cut replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
