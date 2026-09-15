"""Read-only ordinary gravity Ward bridge and entire raw charge report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_massive_matter_graviton_vertex import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-ordinary-gravity-ward-bridge.json"
)
PARENT_SHA = "3018ff11f02ada528818d123caa62a5f552f83c72ff267ed08975b227550c8ab"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_ordinary_gravity_ward_bridge/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError(
            "The frozen matter F1 and conditional background Ward report changed"
        )
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_286_full_matter_F1_conditional_Ward_and_entire_ancestry_rebuilt": PARENT_SHA,
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
        raise ValueError("An ordinary gauge bridge or continuity proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.287.ORIGINAL_ORDINARY_PURE_GRAVITY_ZERO_TRANSFER_WARD_BRIDGE_WITH_WHOLE_TOPOLOGY_SPACELIKE_CONTINUITY_AND_ENTIRE_RAW_CHARGE_LSZ_CANCELLATION",
        "date": "2026-09-15",
        "status": "SCOPED_ORDINARY_GR_CHARGE_AND_POLE_FINITE_CONTINUITY; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/gauge.md",
            "notes/continuity.md",
            "notes/soft.md",
            "notes/matching.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_gauge_bridge": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_continuity_and_raw_charge_normalization": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete on-shell constant-metric gauge-fixing difference is an infrared-integrable scaleless tadpole for the specified linear split and harmonic gauge. Gram-free estimates cover every pure-GR proper-vertex topology and prove spacelike continuity of F1 pole and finite Laurent coefficients. The ordinary zero-transfer proper endpoint and scalar LSZ legs cancel the entire raw residue derivative, including finite7, EulerGamma and log4pi. The remaining nonzero-transfer IR pole has explicit coefficient R(t)/(8pi^2*kappa), which tends to zero but is not a finite physical scattering prescription. Full metric/Newton and curved/four-point local matching, slopes/F2, physical IR/Regge, omitted loops and original-state/bounce obligations remain. Original V/G/B/P8 remain OPEN.",
        "not_established": [
            "Finite physical finite-transfer detector/dressing/IR observable or exact interacting massless LSZ theorem",
            "Full F1 slope, physical F2, H-metric mixing and finite curved Ricci-derivative matching",
            "S285 light pure-curvature polynomial, full physical Newton residue and four-point local/tadpole/double-pole completion",
            "Complete finite crossed amplitude, improved b20, omitted loops or finite Regge contour remainder",
            "Physical EFT cutoff, original quantum state/domain/measure/bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Whole unchanged source, explicit one-loop graph completeness, literal all-D trace and full symmetric-metric GF variation, finite-component independent checks, complete numerator-degree and all-domain Feynman-parameter bounds without Gram division, literal raw triangle sign and entire UV/IR/finite LSZ cancellation. The written analytic and covariance proofs are not FORMALIZED. Original SymPy is retained in native/direct/ordinary/CLI; exact-GCD is FULL-only. No frozen source/raw, finite counterterm or state is edited.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The ordinary gravity Ward bridge report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.287 ordinary gravity Ward bridge and raw charge replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
