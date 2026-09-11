"""Read-only full-source conditional retarded-energy and causal-force replay."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_canonical_affine_decoupling import verify as parent

from . import audit, clock, coefficients, estimates, hamiltonian

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-retarded-energy.json"
PARENT_SHA = "b85c821b59bc41489418c0c190a2cc573cd9852cb19254627c3e2150d1dd6446"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_retarded_energy/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen actual parent and conditional-state chain changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_177_fully_rebuilt": PARENT_SHA,
        "same_fixed_CD_base_parent_and_all_prior_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A full-source conditional response proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.178.COMPLETE_AFFINE_SOURCE_FINITE_TIME_RETARDED_ENERGY_AND_CONDITIONAL_CAUSAL_CUBIC_LIGHT_FORCE",
        "date": "2026-09-11",
        "status": "FULL_SOURCE_FIXED_CD_CONDITIONAL_GAUSSIAN_MEAN_ENERGY_AND_CAUSAL_LIGHT_FORCE; NOT_METRIC_NOISE_FULL_PARENT_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/energy.md",
            "notes/source.md",
            "notes/tame.md",
            "notes/state.md",
            "notes/force.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "physical_retarded_energy_and_full_source_bounds": serialize(
            {
                "hamiltonian": payload(hamiltonian.data()),
                "modes": payload(hamiltonian.modes()),
                "coefficients": payload(coefficients.data()),
                "clock": payload(clock.data()),
                "tame": payload(clock.bounds()),
            }
        ),
        "conditional_mean_energy_and_causal_light_force": serialize(
            payload(estimates.data())
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "On the unchanged CD geometry and the specified compact scalar-clock perturbation class, the complete regular affine source has explicit quadratic Sobolev bounds. The actual sourced Proca temporal constraint and a positive reference Hamiltonian give a finite-time retarded energy estimate at all spatial momenta. With the same conditional Gaussian Cauchy state and canonical measure, this bounds the free coherent mean energy and the full source/contact fixed-metric light force, which is cubically small. Homogeneous closed sources cancel the full force exactly; spatially nonclosed sources are retained. This is not a global on-shell inverse bound, full source-interaction metric tensor, metric/noise response, interacting-parent or physical-cutoff theorem, and V/G/B and original P8 remain open.",
        "not_established": [
            "A metric-response or symmetric-noise norm, variable-geometry response or pointwise total stress bound",
            "Full affine/light/tensor/auxiliary/mixed quantum contributions, physical cutoff or omitted-order control",
            "A self-consistent quantum-corrected background, full interacting state or global-time error bound",
            "Complete vacuum matching/contour/cuts/truncation, finite-gravity IR/Regge remainder or V/G/B closure",
        ],
        "verification_boundary": "Exact constrained Hamiltonian/source identities, full rational/exponential coefficient bounds, explicit fixed-CD jet and Frechet estimates, state/coherent-mean bookkeeping and causal source-contact force norms have written proofs. Independent literal full-source derivatives, nonclosed and homogeneous controls, full-current Fourier evolutions and finite-time resonant tests supplement them. Numerical modes do not replace the all-momentum energy proof or define a physical cutoff. Native, direct science, ordinary and CLI use original SymPy; only full regression uses the audited exact GCD adapter. This is not FORMALIZED or a full quantum-parent certificate.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The full-source conditional retarded-response report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.178.COMPLETE_AFFINE_SOURCE_FINITE_TIME_RETARDED_ENERGY_AND_CONDITIONAL_CAUSAL_CUBIC_LIGHT_FORCE replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
