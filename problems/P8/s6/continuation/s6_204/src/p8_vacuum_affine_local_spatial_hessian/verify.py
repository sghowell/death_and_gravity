"""Read-only fixed local full tracefree spatial Hessian and proper-time bound."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_proca_gaussian import verify as gaussian
from p8_vacuum_affine_reference_spatial_remainder import verify as parent

from . import audit, bounds, helicities, local, proper

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-local-spatial-hessian.json"
PARENT_SHA = "cc0b8a208f63d7a9e170620c4d66c845eea4a362a4c27e252036d7289c975c36"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_local_spatial_hessian/*.py"))
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
        "S6_203_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A fixed local full spatial Hessian gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.204.FIXED_COVARIANT_LOCAL_FULL_TRACEFREE_SPATIAL_HESSIAN_WITH_PROPER_TIME_AND_NORM_CONTROL",
        "date": "2026-09-11",
        "status": "FIXED_COVARIANT_LOCAL_FULL_TRACEFREE_SPATIAL_HESSIAN_WITH_PROPER_TIME_AND_NORM_CONTROL; NOT_QUANTUM_SPATIAL_MATCHING_FULL_MIXED_INVERSE_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/geometry.md",
            "notes/helicities.md",
            "notes/conformal.md",
            "notes/proper-time.md",
            "notes/bounds.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_spatial_polynomial_maps_and_fixed_local_Hessians": serialize(
            {"helicities": payload(helicities.data()), "local": payload(local.data())}
        ),
        "proper_time_conversion_and_actual_CD_derivative_norm": serialize(
            {"proper": payload(proper.data()), "bounds": payload(bounds.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The already-fixed finite local covariant action has an explicit complete Hessian on the actual CD tracefree spatial metric subspace. Smooth polynomial momentum maps retain both tensor, both vector and scalar shear directions without a zero-momentum or light-cone singularity. Literal four-dimensional curvature jets, including noncommuting matrix time derivatives and generic spatial directions, reproduce the full compact Einstein, R-squared and Weyl-squared Hessians and their retained boundary term. The original dimensionally matched finite coefficients are unchanged, and compact Euler variation is removed only in physical dimension after that matching. The exact proper-time density factor a^-1 and all clock/friction terms are retained. The resulting local L2 operator is below5e4 Phi42[Gamma], with two-external-metric canonical coefficient2e-795. This is a fixed local matching target, not proof that the remaining quantum contact and Taylor-integrand sector matches it. No physical state or subtraction prescription is changed, and no full quantum response, reduced mixed inverse, nonlinear background, stability, physical cutoff, remaining parent loops or finite-gravity IR/Regge follows. Original V/G/B and P8 remain open.",
        "not_established": [
            "Fixed covariant matching of the remaining quantum contact/Taylor-integrand sector or its finite/divergent regulator artifacts",
            "A full quantum response or a lapse/shift-constrained, canonically reduced mixed operator",
            "A derivative-compatible inverse, interacting background or stability theorem",
            "Full parent/cutoff/scattering/IR/Regge matching or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact full spatial basis and polynomial maps, complete curvature-Hessian polynomials, original finite coefficients, proper-time measure/clock conversion and actual CD norm constants. Independent literal four-dimensional curvature, all five directions, noncommuting matrix jets, compact boundary terms, generic rotations, wrong-measure controls, self-adjoint pairings and derivative norms supplement written estimates. Native, direct science, ordinary and CLI use original SymPy; only full regression uses the audited exact GCD adapter. The compact geometric variation and continuous norm estimate are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The fixed local full spatial Hessian report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.204.FIXED_COVARIANT_LOCAL_FULL_TRACEFREE_SPATIAL_HESSIAN_WITH_PROPER_TIME_AND_NORM_CONTROL replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
