"""Read-only current-parent conditional prepared homogeneous metric response."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_proca_nonlocal_response import verify as old_response
from p8_vacuum_affine_retarded_energy import verify as parent

from . import audit, bounds, bridge, chart, limits

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-metric-response.json"
PARENT_SHA = "c25ba2599f6297513798ff45291a831aa46be8397dad43fec3d6dc8fdfaea488"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_metric_response/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError(
            "The frozen complete affine conditional-response chain changed"
        )
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    if sha(old_response.REPORT) != bounds.RESPONSE_SHA:
        raise ValueError("The source-pinned ordinary physical response changed")
    old_response.validate_report(
        json.loads(old_response.REPORT.read_text()), old_response.build_report()
    )
    return {
        "S6_178_fully_rebuilt": PARENT_SHA,
        "S6_85_complete_ordinary_response_source_rebuilt": bounds.RESPONSE_SHA,
        "only_old_physical_vector_block_bridged_not_scalar_profile_or_total_clock_response": True,
        "same_full_parent_and_all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A conditional homogeneous metric-response gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.179.SAME_PARENT_CONDITIONAL_HOMOGENEOUS_METRIC_RESPONSE_AND_UNCANCELLED_FULL_CHART_CONTACTS",
        "date": "2026-09-11",
        "status": "ACTUAL_CONDITIONAL_VECTOR_PREPARED_HOMOGENEOUS_METRIC_RESPONSE_BOUNDED_WITH_NONSTATIONARY_CONTACTS; NOT_SPATIAL_NOISE_COUPLED_BACKGROUND_FULL_PARENT_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/bridge.md",
            "notes/homogeneous.md",
            "notes/continuum.md",
            "notes/chart.md",
            "notes/bounds.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "actual_conditional_operator_state_and_homogeneous_source_bridge": serialize(
            {
                "bridge": payload(bridge.data()),
                "homogeneous": payload(bridge.homogeneous()),
            }
        ),
        "complete_prepared_metric_response_and_nonstationary_chart": serialize(
            {
                "chart": payload(chart.data()),
                "chart_envelopes": payload(chart.envelopes()),
                "bounds": payload(bounds.data()),
                "limits": payload(limits.data()),
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged regular affine parent and its specified conditional canonical vector measure/state have the same complete prepared homogeneous physical stress-response operator as ordinary Proca. Literal arbitrary-lapse/scale mode and readout vertices, exact homogeneous source cancellation and the fixed covariant prescription justify transferring the source-pinned all-continuum physical vector norm below5e-795 per kappa. The full analytic chart is differentiated with its nonzero background first variation: readout/volume and second-map contacts give a new clock-current C10-to-C0 norm below1e-784 per kappa. No old scalar counter-profile or background-cancelled clock total is installed. This is not a spatial/noise norm, coupled inverse or self-consistent/full interacting quantum background, and V/G/B and original P8 remain open.",
        "not_established": [
            "Arbitrary spatial or independently prepared initial-state response, a symmetric-noise norm or no-loss coupled inverse",
            "A self-consistent quantum-corrected background, global-time error control or actual total-state tail estimate",
            "Full affine/light/tensor/auxiliary/mixed quantum contributions, interacting measure/state, cutoff or omitted-order control",
            "Complete vacuum matching/contour/cuts/truncation, finite-gravity IR/Regge remainder or V/G/B closure",
        ],
        "verification_boundary": "New exact operator/source/chart/current identities, continuous rational chart envelopes and a source-pinned complete continuum stress bound have written proofs. Every old report is fully rebuilt before certification; only the physically matched vector block is used. Independent literal metric variations, nonlinear chart controls, preparation and derivative-loss tests supplement the argument. Native, direct science, ordinary and CLI use original SymPy; only full regression uses its audited exact GCD adapter. Not FORMALIZED or a full quantum-parent certificate.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The conditional homogeneous metric-response report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.179.SAME_PARENT_CONDITIONAL_HOMOGENEOUS_METRIC_RESPONSE_AND_UNCANCELLED_FULL_CHART_CONTACTS replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
