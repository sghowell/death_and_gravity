"""Read-only state-transport connector and two-real matched soft reference."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_single_residual_soft_dressing import verify as dressing_input
from p8_vacuum_affine_two_real_soft_overlap import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-state-transport-soft-matching.json"
)
PARENT_SHA = "79539dcf68142c2ebb62a0236293be98d2f1b6eece93d49e46cbf091fb8e1994"
DRESSING_SHA = "22e89710b44e3df51ab35cf4ef541f70c42c68ee740a29f261f9f9d16694b59c"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_state_transport_soft_matching/*.py"))
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
            "The frozen state-correct two-real soft-overlap parent changed"
        )
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(dressing_input.REPORT) != DRESSING_SHA:
        raise ValueError(
            "The frozen single-residual leading-soft dressing input changed"
        )
    dressing_input.validate_report(
        json.loads(dressing_input.REPORT.read_text()), dressing_input.build_report()
    )
    return {
        "S6_313_uniform_state_correct_two_real_overlap_rebuilt": PARENT_SHA,
        "S6_309_single_residual_leading_soft_dressing_rebuilt": DRESSING_SHA,
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
        raise ValueError("A state-transport matching proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.314.STATE_TRANSPORT_CONNECTOR_AND_TWO_REAL_MATCHED_SOFT_REFERENCE",
        "date": "2026-09-16",
        "status": "SCOPED_STATE_TRANSPORT_AND_TWO_REAL_MATCHED_SOFT_REFERENCE; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/continuity.md",
            "notes/series.md",
            "notes/matching.md",
            "notes/reference.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_state_continuity": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_signed_series_and_matched_reference": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "Total marked-energy continuity of the full angular conversion and entire signed-series derivative bounds make a necessary state-transport connector finite. Probability-level two-real subtraction and its two-marked Bose dressing then match the complete original0/1/2-real physical D4 tree densities in the stated leading-soft reference. Its relative difference from P0 is below1e-653 and vanishes uniformly at zero threshold. This is not finite hard real-virtual matching, an all-N nonleading inclusive rate, unitarity or V/G/B/P8 closure.",
        "not_established": [
            "Complete finite hard real-virtual and evanescent matching",
            "All-N nonleading real seeds, hard loops or an interacting quantum state",
            "A positive event measure, threshold monotonicity or unitarity",
            "Absolute complex Regge or the original common-parent bounce",
            "Original V/G/B/P8 closure",
        ],
        "verification_boundary": "General complex probability and real-density matching identities; a nuclear-norm and joint radial Holder state-continuity proof; global inverse-Gamma derivative and entire signed-series bounds; exact moment/Bose identities;36 actual recoil-current calibrations; independent series derivatives and connector quadrature; six unchanged original434/47-face decompositions. Finite samples are supplementary, not the uniform proof. Native/direct/ordinary/CLI retain original SymPy; the adapter is FULL-only. All frozen sources and scoped P8(a) are unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The state-transport matched soft-reference report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.314 state-transport matched soft-reference replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
