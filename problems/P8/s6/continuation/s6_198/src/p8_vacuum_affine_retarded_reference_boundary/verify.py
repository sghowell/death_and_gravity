"""Read-only retarded reference boundaries and full finite bulk."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_proca_gaussian import verify as gaussian
from p8_vacuum_affine_retarded_state_remainder import verify as parent

from . import audit, boundary, bounds, ibp, majorants

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-retarded-reference-boundary.json"
)
PARENT_SHA = "d1cbd532c2d04c2875f6c20a7d8cb206ce07815da07e11f84442aa626b32a957"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_retarded_reference_boundary/*.py"))
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
        "S6_197_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A retarded reference boundary gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.198.ACTUAL_SPATIAL_REFERENCE_SIX_STEP_RETARDED_BOUNDARY_EXTRACTION_AND_FULL_CONTINUUM_FINITE_BULK",
        "date": "2026-09-11",
        "status": "ACTUAL_SPATIAL_REFERENCE_SIX_STEP_RETARDED_BOUNDARY_EXTRACTION_AND_FULL_CONTINUUM_FINITE_BULK; NOT_CONTACT_ENDPOINT_COVARIANT_SPATIAL_MATCHING_FULL_INVERSE_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/bounds.md",
            "notes/ibp.md",
            "notes/boundary.md",
            "notes/majorants.md",
            "notes/limit.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "exact_retarded_identity_and_Cauchy_coefficients": serialize(
            {"ibp": payload(ibp.data()), "majorants": payload(majorants.data())}
        ),
        "finite_reference_and_full_remaining_sector": serialize(
            {"bounds": payload(bounds.data()), "boundary": payload(boundary.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The actual constant-alpha W8 reference has an exact six-step retarded integration identity with all initial and equal-time terms accounted for. Complete Cauchy-Leibniz coefficients through sixth order, retaining the full constrained pair amplitude, give a convergent sixth bulk and fifth endpoint uniformly in all external spatial momenta. Their full continuum bound is1e48 in detector L2 and source time-H6 norms, with quantitative two-leg regulator error1e52/K. The corresponding fifth-bulk absolute majorant is not integrable. Adding the already controlled S197 actual-state memory/contact remainder gives a known finite part below2e48 in compatible six-time/one-space source and spatial-H1 detector norms, with canonical display8e-752 and tail8e-748/K. The remaining full reference metric contact and five equal-time endpoints retain their different regulator regions and still require the original covariant spatial matching. Equal-time source jets do not imply spatial locality or a reference Ward identity. No counterterm is retuned, no odd endpoint is discarded, and no full inverse, interacting background, physical cutoff or original V/G/B and P8 closure follows.",
        "not_established": [
            "Original covariant spatial matching of the full reference contact and five equal-time endpoint kernels",
            "Automatic spatial locality, cancellation of odd endpoints or an independent physical reference Ward identity",
            "Full spatial/mixed response inverse, finite-amplitude interacting background, stability or physical cutoff",
            "Full parent state/measure/loops, vacuum cuts/contour/truncation, finite-gravity IR/Regge or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact six-step retarded ODE identities, complete Cauchy-Leibniz coefficients, full-continuum and removed-union radial bounds and canonical factors. Independent literal variable-phase integrals, rational derivative fixtures, lower-boundary, odd-endpoint and nonpolynomial spatial controls supplement written analytic arguments. Native, direct science, ordinary and CLI use original SymPy; only full regression uses the audited exact GCD adapter. The analytic continuity and convergence arguments are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The retarded reference boundary report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.198.ACTUAL_SPATIAL_REFERENCE_SIX_STEP_RETARDED_BOUNDARY_EXTRACTION_AND_FULL_CONTINUUM_FINITE_BULK replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
