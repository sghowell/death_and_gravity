"""Read-only spatial two-momentum current and finite-band response."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_covariant_current import verify as parent
from p8_vacuum_affine_proca_gaussian import verify as gaussian

from . import audit, bounds, hamiltonian, response, vertices

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-spatial-current.json"
PARENT_SHA = "16b91d124384e673e03edb32e2303e5572aeac2fe106593151bfb9b6e8550d8f"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_spatial_current/*.py"))
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
        "S6_194_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A spatial finite-band current response gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.195.ACTUAL_SPATIAL_TWO_MOMENTUM_PROCA_METRIC_VERTICES_CONTACT_AND_UNIFORM_FINITE_BAND_GAUSSIAN_RESPONSE",
        "date": "2026-09-11",
        "status": "ACTUAL_SPATIAL_TWO_MOMENTUM_PROCA_METRIC_VERTICES_CONTACT_AND_UNIFORM_FINITE_BAND_GAUSSIAN_RESPONSE; NOT_RENORMALIZED_SPATIAL_CONTINUUM_INVERSE_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/hamiltonian.md",
            "notes/vertices.md",
            "notes/response.md",
            "notes/bounds.md",
            "notes/limit.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "full_spatial_Hamiltonian_and_vertices": serialize(
            {
                "hamiltonian": payload(hamiltonian.data()),
                "vertices": payload(vertices.data()),
            }
        ),
        "two_momentum_response_and_finite_band": serialize(
            {"response": payload(response.data()), "bounds": payload(bounds.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged actual conditional Proca sector has the full spatially varying unimodular Hamiltonian with its positive temporal-constraint energy. Exact first metric vertices retain two distinct momenta, while the second variation retains the pointwise ordered contact and its full Fourier convolution. A ten-row energy factorization bounds the two-sided vertices uniformly for all momenta and preserves the constrained longitudinal mode. The exact off-diagonal covariance tangent has two distinct propagators and gives the same finite-regulator sign and factor as the Kubo commutator, together with the local metric contact. At computational K1e16 the complete projected current response is bounded by2e104 uniformly in external momentum, giving canonical display8e-696. The contact can survive when external momentum exceeds2K and the propagation overlap is empty. The internal-momentum majorants grow and are not a renormalized spatial continuum bound. Full spatial subtraction/tail matching, response inverse, finite-coupling interacting background, physical cutoff and original V/G/B and P8 remain open.",
        "not_established": [
            "A covariantly renormalized full spatial/mixed continuum response and quantitative infinite-tail matching",
            "A full feedback inverse, finite-coupling interacting quantum background or its stability",
            "Full parent state/measure/loops, physical heavy/cutoff/threshold and omitted-order matching",
            "Complete vacuum cuts/contour/truncation, finite-gravity IR/Regge or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact constrained spatial Hamiltonian and two-momentum feature identities, full local contacts, covariance/Kubo algebra and quantitative finite-band majorants. Independent real-space Legendre, noncommuting lattice Fourier variations and direct coupled Gaussian evolution supplement continuous written proofs. Native, direct science, ordinary and CLI use original SymPy; only full regression uses the audited exact GCD adapter. Analytic energy and continuum-scope arguments are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The spatial finite-band current report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.195.ACTUAL_SPATIAL_TWO_MOMENTUM_PROCA_METRIC_VERTICES_CONTACT_AND_UNIFORM_FINITE_BAND_GAUSSIAN_RESPONSE replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
