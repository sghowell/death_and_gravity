"""Read-only complete fixed covariant homogeneous current and canonical response."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_adiabatic_action import verify as parent
from p8_vacuum_affine_proca_gaussian import verify as gaussian

from . import audit, canonical, current, metrics, tensors

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-covariant-current.json"
PARENT_SHA = "98888c7a7fdcfe4fee142833a4220c62192b2834f736625192c8534201dca278"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_covariant_current/*.py"))
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
        "S6_193_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A complete fixed covariant current response gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.194.ACTUAL_COMPLETE_FIXED_COVARIANT_HOMOGENEOUS_CURRENT_C2_REMAINDER_AND_CANONICAL_RESPONSE",
        "date": "2026-09-11",
        "status": "ACTUAL_COMPLETE_FIXED_COVARIANT_HOMOGENEOUS_CURRENT_C2_REMAINDER_AND_CANONICAL_RESPONSE; NOT_FULL_SPATIAL_MIXED_RESPONSE_INVERSE_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/metrics.md",
            "notes/tensors.md",
            "notes/local.md",
            "notes/current.md",
            "notes/canonical.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_metric_and_local_Euler_bounds": serialize(
            {"metrics": payload(metrics.data()), "tensors": payload(tensors.data())}
        ),
        "full_current_and_canonical_response": serialize(
            {"current": payload(current.data()), "canonical": payload(canonical.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "For the unchanged actual Gaussian Proca state and fixed covariant prescription, the complete homogeneous local current and its first two amplitude derivatives are bounded by1e86. The proof uses full factorial-weighted absolute metric jets, all covariant spatial connections, the exact curvature Euler tensors and both raised metric indices and detector Jacobian. Combining the S193 fixed-prescription identity with the S192 convergent current comparison gives full derivative displays1e87,1e95,2e111 and a finite-amplitude Taylor remainder at most epsilon^2*1e111, with its zero edge retained. The isotropic clock current vanishes by full three-mode rotation invariance. Both canonical field chain factors give a linear homogeneous C12-to-C0 response bound4e-705. No derivative-losing bound is promoted to a same-space contraction or inverse; full spatial/mixed feedback, finite-coupling interacting background, physical cutoff and original V/G/B and P8 remain open.",
        "not_established": [
            "A complete spatial/mixed physical response or general nonlinear Frechet norm beyond the stated amplitude lines",
            "A full feedback inverse, finite-coupling interacting quantum background or its stability",
            "Full parent state/measure/loops, physical heavy/cutoff/threshold and omitted-order matching",
            "Complete vacuum cuts/contour/truncation, finite-gravity IR/Regge or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact rational mixed-jet and tensor majorants, full curvature Euler formulas, local finite-action coefficients, complete-current sums, Taylor integral and canonical normalization. Independent literal spatial metric-jet Euler variations and noncommuting exponential derivatives supplement continuous written proofs. Native, direct science, ordinary and CLI use original SymPy; only full regression uses the audited exact GCD adapter. Analytic estimates and differentiability arguments are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete fixed covariant current report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.194.ACTUAL_COMPLETE_FIXED_COVARIANT_HOMOGENEOUS_CURRENT_C2_REMAINDER_AND_CANONICAL_RESPONSE replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
