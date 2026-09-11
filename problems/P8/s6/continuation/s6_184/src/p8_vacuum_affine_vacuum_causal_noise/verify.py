"""Read-only complete vacuum source, causal mean and scalar-force noise report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_classical_propagator import verify as parent
from p8_vacuum_affine_proca_gaussian import verify as gaussian

from . import audit, energy, noise, source, vacuum

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-vacuum-causal-noise.json"
PARENT_SHA = "d0809a0175afeb864fc015cc59e2b55c6364d0ac8b0982b9118aa115b1268234"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_vacuum_causal_noise/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA or sha(gaussian.REPORT) != GAUSSIAN_SHA:
        raise ValueError("The source-pinned actual affine and Gaussian chain changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_183_fully_rebuilt": PARENT_SHA,
        "S6_176_same_canonical_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "Minkowski_vacuum_already_specified_by_S6_182": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A complete vacuum source response gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.184.COMPLETE_CANONICAL_VACUUM_AFFINE_SOURCE_FINITE_TIME_CAUSAL_MEAN_SCALAR_FORCE_NOISE_AND_CONDITIONAL_GAUSSIAN_DECOUPLING",
        "date": "2026-09-11",
        "status": "COMPLETE_SOURCE_FINITE_TIME_CONDITIONAL_GAUSSIAN_MEAN_AND_SCALAR_FORCE_NOISE_BOUNDS; NOT_FULL_INTERACTING_QUANTUM_MATCHING_METRIC_NOISE_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/tame.md",
            "notes/energy.md",
            "notes/noise.md",
            "notes/vacuum.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_source_and_causal_mean": serialize(
            {"source": payload(source.data()), "mean": payload(energy.data())}
        ),
        "full_scalar_force_noise_and_fixed_vacuum_family": serialize(
            {"noise": payload(noise.data()), "vacuum": payload(vacuum.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete regular canonical affine source admits a uniform complex-jet and Sobolev bound, without discarding its infinitely many harmonics or null-gradient terms. The actual temporal constraint and all-momentum finite-time energy estimate bound its full causal scalar mean force below1e-2378 on the prepared unit-jet vacuum class at the anchor. All three Minkowski Gaussian polarizations give smeared scalar-force standard deviation below2e-1190 in the stated test norm. The full fixed scalar profile and force are below1e-818500 after the explicit same-prescription vacuum constant cancellation. At fixed complete canonical functions and vector mass, source mean and noise variance decay as1/kappa. This is conditional sector control, not an interacting vacuum amplitude, metric-noise or quantum-background bound, physical cutoff or original V/G/B closure.",
        "not_established": [
            "Interacting light, tensor, auxiliary or mixed-loop quantum measures and scattering amplitudes",
            "Metric stress-tensor noise, a stochastic or nonlinear light/gravity solution, or corrected all-time stability",
            "Physical heavy/cutoff, threshold, omitted-order and complete common-parent matching control",
            "Complete vacuum cuts/contour/truncation, finite-gravity IR/Regge or original V/G/B closure",
        ],
        "verification_boundary": "Exact full canonical source and constrained Hamiltonian identities, continuous complex-disc and all-momentum energy bounds, complete physical polarization covariance and full fixed-profile exponent estimates. Independent complex coefficient, source-derivative, constrained-mode, spectral covariance and nonconstant profile fixtures supplement the written proof. Native, direct science, ordinary and CLI use original SymPy; only full regression uses the audited exact GCD adapter. Analytic and Gaussian continuum arguments are not FORMALIZED; no all-frequency physical EFT claim follows.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete vacuum causal/noise report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.184.COMPLETE_CANONICAL_VACUUM_AFFINE_SOURCE_FINITE_TIME_CAUSAL_MEAN_SCALAR_FORCE_NOISE_AND_CONDITIONAL_GAUSSIAN_DECOUPLING replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
