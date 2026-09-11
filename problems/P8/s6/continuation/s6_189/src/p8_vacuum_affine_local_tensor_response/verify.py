"""Read-only fixed finite local full-spatial tensor response Hessian."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_proca_gaussian import verify as gaussian
from p8_vacuum_affine_shear_response import verify as parent

from . import audit, bounds, geometry, operator, prescription

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-local-tensor-response.json"
PARENT_SHA = "017e45a3054f3db472ce50ac88c0bc33328127c4d970a184b18552c0e6789f99"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_local_tensor_response/*.py"))
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
        "S6_188_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("An actual fixed local tensor-response gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.189.ACTUAL_FIXED_PROCA_FINITE_LOCAL_COVARIANT_FULL_SPATIAL_TENSOR_HESSIAN_AND_CANONICAL_H4_BOUND",
        "date": "2026-09-11",
        "status": "ACTUAL_FIXED_FINITE_LOCAL_FULL_SPATIAL_TENSOR_HESSIAN_AND_H4_BOUND; NOT_COMPLETE_DETERMINANT_RESPONSE_INVERSE_FINITE_COUPLING_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/curvature.md",
            "notes/prescription.md",
            "notes/operator.md",
            "notes/bounds.md",
            "notes/contact.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "direct_curvature_and_fixed_prescription": serialize(
            {
                "curvature": payload(geometry.data()),
                "prescription": payload(prescription.data()),
            }
        ),
        "full_spatial_operator_and_Sobolev_bound": serialize(
            {"operator": payload(operator.data()), "bounds": payload(bounds.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The same source-pinned finite local Proca matching density has a complete compact-spacetime TT metric Hessian with the physical canonical normalization. Direct dynamic four-dimensional curvature, the compact Weyl boundary and independent lapse/scale stress variations fix its signs and factors. Its operator retains the A-weighted wave term and full fourth-order Ddag D term with all spatial derivatives and the actual a^3 weighted adjoint. On the unchanged CD slab, its canonical H4-to-L2 norm is below6e-796. This is the Hessian of the specified finite local action, not the entire state/subtraction-dependent determinant response or a small full inverse, finite-coupling error or physical cutoff. Original V/G/B and P8 remain OPEN.",
        "not_established": [
            "A complete state/subtraction-dependent determinant response or uniform integrable UV-subtracted response remainder",
            "All remaining renormalized contacts, a full response inverse, finite-coupling feedback or nonlinear quantum background",
            "A physical cutoff or pole-deletion prescription from the small fixed local operator coefficient",
            "Full interacting parent matching, finite-gravity IR/Regge or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Direct dynamic four-dimensional curvature, exact local prescription and independent stress-sign variations, full weighted Euler operator and continuous fourth-Sobolev majorants. Independent rotated polarization, compact action-variation and full-momentum norm fixtures supplement the written proofs. Native, direct science, ordinary and CLI use original SymPy; only full regression uses the audited exact GCD adapter. Analytic and functional arguments are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The actual fixed local tensor-response report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.189.ACTUAL_FIXED_PROCA_FINITE_LOCAL_COVARIANT_FULL_SPATIAL_TENSOR_HESSIAN_AND_CANONICAL_H4_BOUND replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
