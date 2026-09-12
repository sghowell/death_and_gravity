"""Read-only canonical flat full tensor cut and dispersion benchmark."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_proca_gaussian import verify as gaussian
from p8_vacuum_affine_retarded_reference_boundary import verify as parent

from . import audit, curvature, dispersion, polarizations, projectors

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-flat-tensor-cut.json"
PARENT_SHA = "2aade0843ba6f18fa9bbd6352b66fbf0a95c22d49949f56d90990d15f9d8ebe2"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_flat_tensor_cut/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    for previous, expected in ((parent, PARENT_SHA), (gaussian, GAUSSIAN_SHA)):
        if sha(previous.REPORT) != expected:
            raise ValueError(
                "The source-pinned actual affine Gaussian state chain changed"
            )
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_198_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A canonical flat tensor cut gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.199.ACTUAL_CANONICAL_FLAT_PROCA_FULL_TENSOR_CUT_SUBTRACTED_DISPERSION_AND_CURVATURE_WEIGHT_BENCHMARK",
        "date": "2026-09-11",
        "status": "ACTUAL_CANONICAL_FLAT_PROCA_FULL_TENSOR_CUT_SUBTRACTED_DISPERSION_AND_CURVATURE_WEIGHT_BENCHMARK; NOT_CURVED_SPATIAL_MATCHING_FULL_PARENT_AMPLITUDE_INVERSE_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/dispersion.md",
            "notes/polarizations.md",
            "notes/curvature.md",
            "notes/projection.md",
            "notes/limit.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "full_covariant_stress_and_tensor_spectral_cut": serialize(
            {
                "polarizations": payload(polarizations.data()),
                "projectors": payload(projectors.data()),
            }
        ),
        "controlled_dispersion_and_curvature_benchmark": serialize(
            {
                "dispersion": payload(dispersion.data()),
                "curvature": payload(curvature.data()),
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The canonical mass1000 Proca flat reference vacuum has a full conserved stress cut derived from all nine physical polarization pairs, including the temporal constraint and longitudinal mode. Complete angular projection gives both spin2 and spin0, with densities beta(13s^2+56m^2s+48m^4)/(3840pi^2) and beta(s^2-4m^2s+12m^4)/(384pi^2) for the actual i/4 retarded current. They vanish through the massive threshold and are strictly positive above it. Three subtractions define a convergent dispersive part with exact positive moments, low-energy remainder bounds and controlled infinite spectral tails. Independent four-dimensional curvature Hessians match the ultraviolet tensor weights13/120 and1/72, retaining the Euler term and longitudinal contribution. Both canonical tensor factors multiply the density by4/kappa without providing uniform high-energy smallness. The absorptive cut does not fix finite local polynomials, redo the dimensionally continued prescription, match the curved CD contact/endpoint sector, or establish full parent scattering cuts, a response inverse, interacting background, physical cutoff or finite-gravity IR/Regge control. Original V/G/B and P8 remain open.",
        "not_established": [
            "Original finite local tensor/contact matching or the curved CD reference contact/endpoint response",
            "Full interacting parent vacuum scattering cuts, truncation errors or finite-gravity IR/Regge data",
            "Full response inverse, finite-coupling interacting background, stability or physical cutoff matching",
            "A massless Maxwell replacement, a uniformly small all-energy response, or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact full physical stress pairs, conserved tensor projectors, phase-space/Wick/current factors, spectral moments and ultraviolet curvature weights. Independent covariant potentials, noncollinear boosts, full angular tensors, oscillator current signs, spectral integrals and linearized four-dimensional curvature Hessians supplement written arguments. Native, direct science, ordinary and CLI use original SymPy; only full regression uses the audited exact GCD adapter. The analytic dispersion and covariance arguments are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The canonical flat tensor cut report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.199.ACTUAL_CANONICAL_FLAT_PROCA_FULL_TENSOR_CUT_SUBTRACTED_DISPERSION_AND_CURVATURE_WEIGHT_BENCHMARK replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
