"""Read-only complete mixed heavy-gravity amplitude and known forward bound."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_whole_quartic_gravity_sector import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-whole-mixed-heavy-gravity-sector.json"
)
PARENT_SHA = "ad67bca9916c2ec0a094f28b1783bd05e070c7fddd345ffead3d6c61bc7a64d8"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_whole_mixed_heavy_gravity_sector/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen whole quartic-gravity parent report changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_293_whole_quartic_gravity_and_known_bound_entire_ancestry_rebuilt": PARENT_SHA,
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
        raise ValueError("A complete mixed-amplitude/forward proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.294.COMPLETE_MINIMAL_MIXED_HEAVY_GRAVITY_WHOLE_D_OFFSHELL_AMPLITUDE_ALL_TWO_BODY_CUTS_FINITE_FORWARD_AND_EXPLICIT_KNOWN_BOUND",
        "date": "2026-09-15",
        "status": "SCOPED_COMPLETE_KNOWN_MIXED_GRAVITY_COEFFICIENT_AND_EXPLICIT_BOUND; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/graphs.md",
            "notes/cuts.md",
            "notes/forward.md",
            "notes/matching.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_complete_mixed_graphs": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_mixed_cuts_and_finite_forward_coefficient": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The entire selected minimal g^2/kappa graph sector includes twelve boxes, full offshell heavy endpoints, selfenergy, four external residues, the required two metric-cubic contact bubble and complete inherited matter-metric endpoint. Whole-D cuts and the exact infrared-subtracted finite amplitude have an explicit original known forward bound below10^-1001. Independent R H/heavy-residue and higher matching, other coupling/all-loop terms, physical IR/Regge and original V/G/B/P8 remain open.",
        "not_established": [
            "Independent finite R H, additional heavy-residue or higher-EFT matching",
            "Positive physical spectral measure or isolated forward positivity verdict",
            "Stable heavy atom, physical width resummation or all-channel dressed IR",
            "Fixed-transfer Regge, other coupling or all-loop controlled remainder",
            "Original state/domain/measure/bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Complete labelled graph inventory, general-D and literal component contractions, whole physical tree cuts, full covariant R H dictionary, exact finite box subtraction and explicit positive-domain majorants support this scoped known loop result. Independent whole-D finite quadratures are calibrations, not all-domain formal proofs or rigorous quadrature-error certificates. The inherited analytic soft division remains conditional. Native/direct/ordinary/CLI use original SymPy; guarded exact-GCD is FULL-only. Frozen sources/raw and original closure status are unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete mixed-gravity/forward report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.294 complete mixed heavy-gravity amplitude and known forward bound replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
