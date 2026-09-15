"""Read-only whole original graviton normal-cut sheet and pole-running report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_massive_elastic_proca_infrared import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-massive-graviton-threshold-sheet.json"
)
PARENT_SHA = "7cbf5dae85909e1e562ddce6dba70ea0646c6afecbcd552d3b48b9fda22c865b"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_massive_graviton_threshold_sheet/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen original elastic and Proca cut changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_281_complete_elastic_Proca_cut_and_entire_ancestry_rebuilt": PARENT_SHA,
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
        raise ValueError("A massive graviton normal-cut proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.282.ORIGINAL_COMPLETE_MASSIVE_GRAVITON_NORMAL_CUT_INVARIANT_SHEET_NONUNIFORM_THRESHOLD_AND_REQUIRED_SUBTRACTION_POLE_RUNNING",
        "date": "2026-09-15",
        "status": "SCOPED_COMPLETE_NORMAL_GRAVITON_CUT_SHEET_AND_POLE_RUNNING; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/angular.md",
            "notes/sheets.md",
            "notes/threshold.md",
            "notes/subtraction.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_closed_angular_cut": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_invariant_sheet_threshold_and_pole_running": serialize(
            {
                "packets": {name: payload(packets[name]) for name in names[2:]},
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "Exact original four-helicity normalization and complete two-axis rational integration give the whole massive graviton normal cut. Its invariant representation uses the same analytic massive pair kernel as S278 and removes the apparent external-threshold pole by a uniform integral. The physical normal-cut branch continues explicitly, including crossed massive boundaries; it is not equated to the full imaginary part in overlapping-cut regions. Complete forward and strict-angle massless endpoint expansions are nonuniform. A convergent pole-factor-preserving channel lift requires cprime(L)=-rho(L)/pi when only the ordinary upper shell changes; its initial physical pole anchor remains undetermined. Dropping that anchor changes the crossing-even b20 coefficient. Full real crossed matching, transfer subtraction, finite local/detector/Regge errors and original quantum-state/bounce completion remain open. Original V/G/B/P8 remain OPEN.",
        "not_established": [
            "Exact interacting LSZ vacuum, original finite massless-pole anchor or all-orders factorization",
            "Complete crossed or double-spectral amplitude, transfer low-cut subtraction or forward b20",
            "Unknown real local matching, omitted loops, finite detector errors or finite Regge remainder",
            "A physical EFT cutoff, unique full off-background dimensional continuation or UV completion",
            "Original quantum state/domain/measure/bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Whole unchanged source and literal four-helicity parent normalization; exact full two-axis division, parameter primitives, invariant kernel identities and uniform threshold remainder. Written complex analysis establishes the declared normal-sheet continuation and endpoint Cauchy construction, not a full amplitude or unknown pole normalization. Written proofs are not FORMALIZED. Numerical diagnostics do not prove continuum statements. Native/direct/ordinary/CLI retain original SymPy; the exact-GCD adapter is restricted to a captured complete FULL regression. Frozen science and all historical qualifications remain unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete graviton threshold sheet report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.282 complete original graviton normal-cut sheet and pole-running replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
