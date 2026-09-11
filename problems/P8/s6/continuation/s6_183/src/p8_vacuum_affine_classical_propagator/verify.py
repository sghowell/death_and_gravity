"""Read-only actual classical propagator and first formal response certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_metric_response import verify as metric
from p8_vacuum_affine_quantum_retuning import verify as parent

from . import audit, hamiltonian, jet_bounds, perturbation, real_bounds

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-classical-propagator.json"
PARENT_SHA = "f754f6d62795d5a1c0ade4c576395f1587fc35f9c9895bc8b7b15b2b2dec2e16"
METRIC_SHA = "17a635f2d77b93b8a815cc3e89eae976c06e80d1880c1df02b475d41aec2ad78"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_classical_propagator/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA or sha(metric.REPORT) != METRIC_SHA:
        raise ValueError("The actual retuned mean or physical response chain changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_182_fully_rebuilt": PARENT_SHA,
        "S6_179_same_physical_response_transitively_rebuilt": METRIC_SHA,
        "actual_S6_180_classical_reference_rederived": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A classical propagator or formal response gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.183.ACTUAL_CLASSICAL_CAUSAL_HAMILTONIAN_C0_AND_C10_INVERSE_AND_FIRST_FORMAL_STATIONARY_GAUSSIAN_RESPONSE_BOUND",
        "date": "2026-09-11",
        "status": "QUANTIFIED_CLASSICAL_PREPARED_PROPAGATOR_AND_FIRST_FORMAL_RESPONSE; NOT_FINITE_COUPLING_ERROR_FULL_INVERSE_STABILITY_QUANTUM_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/Hamiltonian.md",
            "notes/real-norm.md",
            "notes/jets.md",
            "notes/first-order.md",
            "notes/causality.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "actual_Hamiltonian_and_continuous_real_norm": serialize(
            {
                "Hamiltonian": payload(hamiltonian.data()),
                "real_norm": payload(real_bounds.data()),
            }
        ),
        "classical_jet_norm_and_actual_first_formal_response": serialize(
            {
                "C10": payload(jet_bounds.data()),
                "formal_response": payload(perturbation.data()),
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged actual classical reference admits a prepared physical C0 inverse below53000 and C10 inverse below1e96 on the stated slab, with the correct force map and original M1 charge. Uniform complex discs and a scaled derivative recurrence give a continuous ten-jet estimate. The actual new stationary Gaussian-plus-fixed-profile response then gives a unique first formal metric coefficient below4e-669 per unit physical C10 force. Future extension preserves causality without imposing right-endpoint zero data. This is not a finite-coupling error, small full inverse, stability theorem, full quantum background, physical cutoff or original V/G/B closure.",
        "not_established": [
            "An exact solution branch differentiable in the formal loop marker or finite-coupling truncation error",
            "A no-loss quantum endomorphism, C0 feedback contraction or small/stable full retained inverse",
            "Full interacting quantum, nonlinear, spatial/noise, global-time or physical cutoff control",
            "Complete vacuum cuts/contour/truncation, finite-gravity IR/Regge or original V/G/B closure",
        ],
        "verification_boundary": "Exact source-normalized Legendre, weighted Hamiltonian, original Euler and matter-charge identities; continuous real polynomial and complex-disc Cauchy bounds; scaled Leibniz recurrence; actual source-pinned stationary response and norm composition. Independent rational-time, high-precision complex-disc, derivative and ODE fixtures supplement the written proof. Native, direct science, ordinary and CLI use unmodified SymPy; only full regression uses the audited exact GCD adapter. The analytic and formal-parameter arguments are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The classical propagator report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.183.ACTUAL_CLASSICAL_CAUSAL_HAMILTONIAN_C0_AND_C10_INVERSE_AND_FIRST_FORMAL_STATIONARY_GAUSSIAN_RESPONSE_BOUND replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
