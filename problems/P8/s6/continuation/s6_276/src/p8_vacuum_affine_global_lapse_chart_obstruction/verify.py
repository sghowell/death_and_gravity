"""Read-only full-source joint-atlas and connected lapse-chart obstruction report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_canonical_boundary_corrected_hybrid import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-global-lapse-chart-obstruction.json"
)
PARENT_SHA = "4951c18e25f7e51811592467e67e060231f9de7fd09f90722d08b96658fcc537"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_global_lapse_chart_obstruction/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen canonical-boundary corrected report changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_275_corrected_finite_hybrid_and_full_ancestry_rebuilt": PARENT_SHA,
        "S275_corrected_finite_scope_not_refuted": True,
        "S269_to_S274_physical_qualifications_unchanged": True,
        "all_original_primitive_and_matching_frontiers_retained": True,
        "no_new_state_or_global_quantum_dynamics_claim": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A full-source lapse-chart proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.276.COMPLETE_FULL_SOURCE_JOINT_AUXILIARY_ATLAS_AND_CONSTRAINED_FINITE_PHASE_OBSTRUCTION_TO_GLOBAL_SINGLE_LAPSE_CHART",
        "date": "2026-09-15",
        "status": "SCOPED_JOINT_ATLAS_AND_CONNECTED_LAPSE_CHART_OBSTRUCTION; NOT_QUANTUM_OR_ORIGINAL_V_G_B_P8_NO_GO",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/auxiliary.md",
            "notes/fold.md",
            "notes/fields.md",
            "notes/branch.md",
            "notes/validation.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_joint_auxiliary_atlas_and_exact_source_fold": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_constraint_realization_and_connected_sheet_obstruction": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The original joint trace/temporal stationary system is invertible for every positive R; its sequential Gamma zero is not a full stationary singularity. The complete original bounce-source jets, retaining actual fixed quantum profiles and the nonzero primitive derivative, separately produce a genuine lapse fold with N=R=Gamma=1, T=0, C_s=3 and C_NN>28. Exact finite canonical axis-wave inputs preserve all spatial momentum constraints, residual translations and original means, and connect the original local lapse sheet to the fold. Therefore no globally C1 lapse map extends that sheet to every such finite input. The full branch envelope and finite primary/secondary rank loss are explicit. This is an obstruction to this particular global single-sheet extension, not a quantum or original-P8 no-go theorem. The corrected S275 finite hybrid and all earlier physical qualifications remain unchanged; original V/G/B/P8 remain OPEN.",
        "not_established": [
            "A global time solution reaching the fold, Lorentzian singularity or clock-time turning point",
            "A no-go theorem for all constrained quantizations, UV completions or original P8",
            "Automatic new gauge symmetry, a rank-changing quantum measure or unitary branch gluing",
            "A new quantum state, self-adjoint domain, regulator removal or physically matched cutoff",
            "Original UV matching, omitted-loop/Regge control, global completion or formal proof verification",
        ],
        "verification_boundary": "Exact complete-source and full-metric identities, actual fixed-profile enclosures and independent diagnostic fixtures support written analytic proofs. Native/direct/ordinary/CLI require original SymPy; only a separately captured full regression may use the audited exact-GCD adapter. Frozen history is never rewritten. The obstruction concerns a globally smooth continuation of the established local lapse sheet, not the validated finite regulated model or original P8 closure.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete full-source lapse-chart report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.276 full-source joint atlas and connected lapse-chart obstruction replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
