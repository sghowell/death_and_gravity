"""Read-only complete homogeneous flat quadratic one-loop stress."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_flat_local_energy import verify as parent

from . import audit, dimensional, pressure, projector

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-flat-dirac-pressure.json"
PARENT_SHA = "f0d35fab196ee3bc72ebd506d8949ab92c6d6c5d48e6ccb7f5fcaa492c5f4ded"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_flat_dirac_pressure/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen absolute local-energy report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_168_fully_rebuilt": PARENT_SHA,
        "same_exact_free_flat_Hadamard_states_and_mass_reference": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A homogeneous flat stress proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.169.COMPLETE_HOMOGENEOUS_FLAT_DIRAC_ONE_LOOP_STRESS_WITH_DIMENSIONAL_CURVATURE_IMPROVEMENT",
        "date": "2026-09-10",
        "status": "BOUNDED_COMPLETE_HOMOGENEOUS_FLAT_TENSOR_OF_SPECIFIED_QUADRATIC_ONE_LOOP_MS_MR_CONTRIBUTION; NOT_CURVED_INTERACTING_PARENT_HIGHER_LOOPS_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/projector.md",
            "notes/dimensional.md",
            "notes/pressure.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "projector_and_dimensional_pressure": serialize(
            {
                "projector": payload(projector.data()),
                "dimensional": payload(dimensional.data()),
            }
        ),
        "absolute_pressure_and_stress": serialize(payload(pressure.data())),
        "observable_and_source_context": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "Both physical pressure components are retained in the four exact rotating frames. The complete half-line state and projector remainder is below1e410. Full dimensional isotropic and angular factors restore the local MS pressure, including its curvature improvement; the derivative finite part is below1e595. With the same entire potential and fixed mass reference, |P|<1e789 at every real time. The unchanged energy bound controls every component of the specified homogeneous flat tensor in its fixed rest frame. This is not the curved or full interacting parent tensor, a bounce-density relative error, or original P8 closure.",
        "not_established": [
            "Curved Hadamard state, scalar source bound or physical Newton-reference dictionary",
            "Full canonical-field/heavy/gauge parent stress, physical cutoff or higher-loop remainders",
            "Quantum target matching, controlled parent background/response, global V/G/B or original P8 closure",
            "A uniform bound on every boosted tensor component or formalization of the analytic proofs",
        ],
        "verification_boundary": "Exact Pauli/projector and pressure UV algebra, full dimensional gamma finite parts, explicit local improvement and complete radial integrals. Independent full-function frame projectors, finite-angle and 600-digit tiny-angle tests, nonzero-regulator quadratures, full profile derivatives and finite-improvement controls supplement the written estimates. Numerical tests are not validated integration. Every own payload uses the unchanged exact serializer. Native, direct science, ordinary and CLI use original SymPy; only full regression uses its separately audited exact GCD adapter.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The homogeneous flat Dirac stress report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.169.COMPLETE_HOMOGENEOUS_FLAT_DIRAC_ONE_LOOP_STRESS_WITH_DIMENSIONAL_CURVATURE_IMPROVEMENT replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
