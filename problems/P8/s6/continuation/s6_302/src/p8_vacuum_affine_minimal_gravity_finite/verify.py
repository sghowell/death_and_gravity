"""Read-only finite minimal-gravity reference and compact-bound certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_massive_gravity_pole_completion import verify as pole_input
from p8_vacuum_affine_radiative_angular_finite import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-minimal-gravity-finite.json"
PARENT_SHA = "5e18135d1a49c2846047f2a71c3de00f7c82b84a131d8a6f924890308f4ccd99"
POLE_SHA = "658ed46e1d638540bc34c2e2259c2ed3a05ec2fe539005904e2f8342d2f4a543"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_minimal_gravity_finite/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen finite radiative angular parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(pole_input.REPORT) != POLE_SHA:
        raise ValueError("The frozen physical gravity-pole completion changed")
    pole_input.validate_report(
        json.loads(pole_input.REPORT.read_text()), pole_input.build_report()
    )
    return {
        "S6_301_finite_radiative_angular_conversion_rebuilt": PARENT_SHA,
        "S6_288_whole_D_physical_gravity_pole_completion_rebuilt": POLE_SHA,
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
        raise ValueError("A finite minimal-gravity proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.302.COMPLETE_KNOWN_MINIMAL_GRAVITY_FINITE_REFERENCE_THREE_SYMBOLIC_MATCHING_ANCHORS_AND_COMPACT_BOUND",
        "date": "2026-09-15",
        "status": "SCOPED_KNOWN_MINIMAL_GRAVITY_FINITE_REFERENCE_AND_COMPACT_BOUND; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/masters.md",
            "notes/bounds.md",
            "notes/assembly.md",
            "notes/matching.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_finite_master_dictionary": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_finite_gravity_assembly_and_compact_bound": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The known minimal Einstein/massive-scalar one-loop finite reference is completely assembled from unchanged full-D inputs, with all evanescent, Gram, LSZ and soft-division finite terms. Its known compact amplitude is bounded by1e9*(1+abs(log(delta)))/(kappa^2*delta^2). At delta1e-204 its original-parameter known rate interference is below1e-579. Three finite matching coordinates remain unassigned; complete hard and radiative amplitudes, Regge and original V/G/B/P8 remain OPEN.",
        "not_established": [
            "Values or bounds of alpha,beta,delta_kappa",
            "The complete current-action or detector-inclusive amplitude",
            "Uniform exact-forward, high-energy Regge and omitted-loop bounds",
            "Original quantum state/domain/measure/common-parent bounce",
            "Original V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact full-D Laurent, tagged UV/IR, crossing and termwise majorant identities accompany the written Feynman-contour and endpoint bounds. Independent high-precision physical quadratures and direct small-epsilon assemblies calibrate formulas, not rigorous quadrature errors. The bound is on the specified known representative, not on its three unknown matching pieces. Native/direct/ordinary/CLI use original SymPy; the adapter is FULL-only. Frozen ancestors and original closure statuses are unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The finite minimal-gravity report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.302 finite minimal-gravity reference replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
