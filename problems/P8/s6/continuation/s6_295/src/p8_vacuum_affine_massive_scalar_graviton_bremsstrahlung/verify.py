"""Read-only complete scalar bremsstrahlung and physical real soft-error bound."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_whole_mixed_heavy_gravity_sector import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT
    / "certificates/polynomial-vacuum-affine-massive-scalar-graviton-bremsstrahlung.json"
)
PARENT_SHA = "365fbe4a68a52fd00725a6e43341a1c863638289bea1d7b12955067bdc4a22a0"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(
            ROOT.glob(
                "src/p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung/*.py"
            )
        )
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError(
            "The frozen complete mixed heavy-gravity parent report changed"
        )
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_294_complete_mixed_gravity_amplitude_and_bound_entire_ancestry_rebuilt": PARENT_SHA,
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
        raise ValueError("A full scalar bremsstrahlung/real-rate proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.295.COMPLETE26_GRAPH_MASSIVE_SCALAR_BREMSSTRAHLUNG_EXACT_RECOIL_PHASE_SPACE_AND_EXPLICIT_PHYSICAL_REAL_RATE_ERROR",
        "date": "2026-09-15",
        "status": "SCOPED_COMPLETE_SELECTED_REAL_EMISSION_AND_PHYSICAL_SOFT_ERROR_BOUND; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/tree.md",
            "notes/recoil.md",
            "notes/rate.md",
            "notes/matching.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_complete_bremsstrahlung": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_physical_recoil_and_real_rate_error": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The full selected26-graph scalar bremsstrahlung amplitude is Ward consistent and its leading heavy remainder is pure gauge. The exact on-shell recoil/phase-space comparison gives a lower-cutoff-uniform finite real-rate error, bounded by[42248192*resolution+54450000000*resolution^2]/(2pi^2*kappa), below10^-792 at original parameters on the explicit compact domain. Individual real rates still diverge; virtual pairing, pure-gravity radiation, higher matching, full IR/Regge and original V/G/B/P8 remain open.",
        "not_established": [
            "Finite individual Fock real rates or virtual/dressed inclusive observable",
            "Equality to S278 dimensional analytic soft division",
            "Gravity-exchange radiation, hard-loop or higher-EFT matching control",
            "All-energy fixed-transfer Regge or uniform forward gravity limit",
            "Original state/domain/measure/bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "The proof combines a full labelled tree inventory, literal D4/5/6 and general algebraic Ward/TT identities, the exact tuned Born identity, an on-shell recoil boost and positive-domain analytic majorants. Independent finite rate quadratures calibrate the implementation but do not certify integration errors or the all-domain inequality. Both real rates individually retain an IR logarithm. Native/direct/ordinary/CLI use original SymPy; guarded exact-GCD is FULL-only. No frozen source/raw, original matching value or closure status is altered.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete bremsstrahlung/real-rate report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.295 complete scalar bremsstrahlung and physical real-rate error replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
