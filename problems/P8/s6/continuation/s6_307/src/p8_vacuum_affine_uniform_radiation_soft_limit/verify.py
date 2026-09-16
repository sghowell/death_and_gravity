"""Read-only complete single-real uniform soft-threshold certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_minimal_gravity_radiation import verify as radiation_input
from p8_vacuum_affine_mixed_gravity_physical_rate import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-uniform-radiation-soft-limit.json"
)
PARENT_SHA = "dd1e86a53242797536676d21c87acc87965f603ca488d168a7bab5c510dd97f6"
RADIATION_SHA = "a34d5160c55891304eda5f83bf90970854774e8aae327d9375b89e5fa797063e"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_uniform_radiation_soft_limit/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen complete physical mixed-sector parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(radiation_input.REPORT) != RADIATION_SHA:
        raise ValueError("The frozen complete47-tree radiation source changed")
    radiation_input.validate_report(
        json.loads(radiation_input.REPORT.read_text()), radiation_input.build_report()
    )
    return {
        "S6_306_complete_physical_mixed_gravity_parent_rebuilt": PARENT_SHA,
        "S6_304_complete47_tree_single_real_radiation_rebuilt": RADIATION_SHA,
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
        raise ValueError("A uniform radiation soft-limit proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.307.COMPLETE_SINGLE_REAL_UNIFORM_SOFT_LIMIT_AND_RELATIVE_LEADING_REFERENCE_BOUND",
        "date": "2026-09-16",
        "status": "SCOPED_COMPLETE_SINGLE_REAL_UNIFORM_SOFT_LIMIT; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/currents.md",
            "notes/estimates.md",
            "notes/integration.md",
            "notes/softlimit.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_paired_currents": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_uniform_amplitudes_and_soft_threshold": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "Exact pairing of incoming and outgoing external radiation currents, with every off-shell shift retained, improves the complete47-tree low-energy remainder to2e8/sqrt(kappa) independent of hard angle. The full one-graviton real-minus-leading-soft total variation is below(2e14*x+2e37*x^2)/kappa uniformly at all nonforward angles and tends tozero uniformly asx->0. Against the unexpanded elastic known leading-soft factor its relative size is below4e-768 at original parameters and also vanishes uniformly at zero threshold. This is not an all-N nonleading or full-rate theorem. Unknown hard matching and original V/G/B/P8 remain OPEN.",
        "not_established": [
            "Complete nonleading amplitudes and calorimetric recoil uniformly in unbounded graviton multiplicity",
            "Unknown finite physical hard and curvature matching or the complete all-order detector rate",
            "Finite forward cross section or exchanged regulator and resolution limits",
            "Complex Coulomb phase and high-energy Regge estimates or a UV quantum construction",
            "Original V/G/B/P8 closure",
        ],
        "verification_boundary": "General exact stress-shift and paired-current identities, explicit massive recoil gaps and all-diagram bounds, two-region uniform energy integration, and an unexpanded leading-soft reference comparison. Independent exact diagrams and simultaneous forward-soft samples calibrate the estimates but are not the uniform proof. Native/direct/ordinary/CLI use original SymPy; the adapter is FULL-only. Frozen ancestors, scoped P8(a), all matching coordinates and original statuses are unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The uniform radiation soft-limit report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.307 uniform single-real soft-limit replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
