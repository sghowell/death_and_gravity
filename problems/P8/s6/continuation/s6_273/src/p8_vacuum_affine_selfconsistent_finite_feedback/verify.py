"""Read-only complete finite self-consistent hybrid feedback certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_finite_volume_turnaround import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-selfconsistent-finite-feedback.json"
)
PARENT_SHA = "a2a032e79128e374f60ec20a410ff7704314bccad2be44ddbcc8bff2d2d67efd"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_selfconsistent_finite_feedback/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen finite volume-turnaround report changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_272_complete_volume_source_and_ancestry_rebuilt": PARENT_SHA,
        "original_source_state_parameters_and_two_cutoffs_unchanged": True,
        "S261_refutation_and_S265_integrability_boundary_unchanged": True,
        "all_original_primitive_and_matching_frontiers_retained": True,
        "new_model_is_classical_homogeneous_quantum_nonzero_hybrid": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A complete finite hybrid feedback gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.273.COMPLETE_SELFCONSISTENT_FINITE_CLASSICAL_HOMOGENEOUS_QUANTUM_NONZERO_MODE_FEEDBACK_AND_PHYSICAL_VOLUME_TURNAROUND_WITH_UNCHANGED_FULL_SOURCE",
        "date": "2026-09-14",
        "status": "FINITE_SELFCONSISTENT_HYBRID_TURNAROUND; NOT_HOMOGENEOUS_QUANTIZATION_OR_ORIGINAL_V_G_B_P8",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/homogeneous.md",
            "notes/symmetry.md",
            "notes/geometry.md",
            "notes/quantum.md",
            "notes/dynamics.md",
            "notes/validation.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_source_canonical_force_and_actual_reference_symmetry": serialize(
            {name: payload(packets[name]) for name in names[:3]}
        ),
        "whole_full_spatial_average_operator_bounds_and_coupled_turnaround": serialize(
            {name: payload(packets[name]) for name in names[3:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The separately named finite hybrid model couples classical homogeneous density variables to the unchanged48-pair quantum nonzero-mode seed and full original nonlinear source. On |u|<=1e-180, all500 full source derivatives, complete averaged spatial/auxiliary remainders, actual cubic covariance and both operator orderings give a strict coupled path contraction. The unique trajectory stays within1e-390 of the reference homogeneous data. All five homogeneous traceless-shape and three spatial-vector canonical pairs obey their zero equations by symmetry; the heavy scalar and full scalar-phase force remain. The physical-volume relative error is<1e-380, both endpoints exceed the center by>5e-360, and every minimum lies within1e-188. No unique or strict minimum, homogeneous quantization, unlocalized dynamics, original interacting mean, mode/volume limit, matching, omitted-loop/Regge/UV control or nonlinear global completion is inferred. Original V/G/B/P8 remain OPEN.",
        "not_established": [
            "Quantization of homogeneous variables or the original fully quantum mean",
            "A unique or strict positive-acceleration volume minimum",
            "Small state displacement or long-time coherent core leakage",
            "Original unregularized dynamics or mode/volume/cutoff removal",
            "Physical matching, omitted-loop/Regge/UV or global completeness",
            "Formal verification of the written analytic/operator proofs",
        ],
        "verification_boundary": "Exact full source and canonical identities, all finite group actions, outward rational budgets and independent diagnostics support written proofs of the explicitly named finite hybrid model. The entire unchanged ancestry is rebuilt. Native/direct/ordinary/CLI use original SymPy; only a captured complete regression may use the audited exact-GCD adapter. Successful historical replay does not reverse either archived physical erratum or close any original frontier.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete finite hybrid feedback report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.273 complete finite self-consistent HYBRID feedback replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
