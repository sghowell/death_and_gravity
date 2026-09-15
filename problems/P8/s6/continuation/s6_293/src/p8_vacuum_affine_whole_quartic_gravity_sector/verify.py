"""Read-only complete quartic-gravity coefficient and explicit regulator bound."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_heavy_resonance_soft_pole_pairing import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-whole-quartic-gravity-sector.json"
)
PARENT_SHA = "4a9a38e52609567bca2b22c1d1384337e58233e0b30f305d3f67854418c38a8f"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_whole_quartic_gravity_sector/*.py"))
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
            "The frozen minimal heavy-transition/formal-pair parent report changed"
        )
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_292_minimal_heavy_transition_and_formal_soft_pair_entire_ancestry_rebuilt": PARENT_SHA,
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
        raise ValueError("A complete quartic/finite-forward proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.293.COMPLETE_QUARTIC_GRAVITY_WHOLE_D_SCALAR_CUT_CLOSED_FINITE_FORWARD_COEFFICIENT_AND_EXPLICIT_UNIFORM_REGULATOR_BOUND",
        "date": "2026-09-15",
        "status": "SCOPED_COMPLETE_KNOWN_QUARTIC_GRAVITY_COEFFICIENT_AND_EXPLICIT_BOUND; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/proper.md",
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
        "whole_original_source_and_complete_quartic_graphs": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_scalar_cut_and_finite_forward_coefficient": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete minimal one-loop C/kappa graph sum has a full-D massive scalar cut equal to the entire GR tree interference. Its inherited formal S278 soft-divided forward coefficient has a closed Catalan/log expression, an explicit original-domain regulator remainder below8epsilon*absC/kappa, and a known magnitude below10^-1005. Constant quartic/curvature anchors have zero b20 but higher matching, all other coupling sectors, physical IR/Regge and original V/G/B/P8 remain open.",
        "not_established": [
            "Complete higher-EFT or other-coupling full-source matching",
            "Positive physical spectral measure or isolated b20 positivity verdict",
            "Exact heavy width, all-channel physical Coulomb/detector/dressed IR",
            "Physical soft-factor unitarity, fixed-transfer Regge or all-loop errors",
            "Original state/domain/measure/bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Complete graph inventory, full D tensor and OS reduction, entire scalar-cut/tree sewing, exact convergent center integrals and written uniform majorants support the scoped result. Whole-D numerical checks are independent tests, not all-domain formal proofs. The original S278 analytic soft-factor convention is retained conditionally; no physical detector or finite matching is selected. Native/direct/ordinary/CLI use original SymPy, with guarded exact-GCD FULL-only. No frozen source/raw or original closure status changes.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete quartic/finite-forward report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.293 complete quartic-gravity coefficient and explicit regulator bound replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
