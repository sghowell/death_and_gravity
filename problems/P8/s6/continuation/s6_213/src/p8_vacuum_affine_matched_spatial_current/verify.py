"""Read-only complete matched tracefree spatial Gaussian current."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_dimensional_spatial_symbol import verify as parent
from p8_vacuum_affine_proca_gaussian import verify as gaussian

from . import assembly, audit, dimension, lift, regulator

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-matched-spatial-current.json"
PARENT_SHA = "151f79113368b3c3c704fc9a7479696b115e24ea5d26bacc0e084083f7b70138"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_matched_spatial_current/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {key: value for key, value in d.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    for previous, expected in ((parent, PARENT_SHA), (gaussian, GAUSSIAN_SHA)):
        if sha(previous.REPORT) != expected:
            raise ValueError(
                "The source-pinned actual affine Gaussian state chain changed"
            )
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_212_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A complete matched tracefree spatial current gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.213.COMPLETE_MATCHED_TRACEFREE_SPATIAL_GAUSSIAN_CURRENT_AND_ANCHORED_UNIFORM_REGULATOR_LIMIT",
        "date": "2026-09-11",
        "status": "COMPLETE_MATCHED_TRACEFREE_SPATIAL_GAUSSIAN_CURRENT_AND_ANCHORED_UNIFORM_REGULATOR_LIMIT; NOT_FULL_SCALAR_MIXED_INVERSE_FINITE_COUPLING_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/dimension.md",
            "notes/bounds.md",
            "notes/anchor.md",
            "notes/lift.md",
            "notes/regulator.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "full_dimensional_limit_and_homogeneous_Hilbert_lift": serialize(
            {"dimension": payload(dimension.data()), "lift": payload(lift.data())}
        ),
        "complete_matched_current_and_original_anchored_regulator": serialize(
            {
                "assembly": payload(assembly.data()),
                "regulator": payload(regulator.data()),
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete tracefree spatial Gaussian retarded current is assembled in the unchanged original prescription. All four actual dimension-dependent WKB coefficients give explicit nonzero frequency bounds on|d-3|<=1/4. Full general-dimensional angular geometry, unexpanded low band, near complement and far Cauchy remainder supply a uniform integrable all-transfer majorant; the subtracted comparison dimension limit equals the physical remainder. The exact common prepared covariance tangent and full contact identify the fixed homogeneous anchor. Original finite matching gives Jren=H0+Known(P)-Known(0)+Ffinite(P), retaining every state/time/endpoint contribution and all finite evanescent terms without double counting. Tensor duality, complex Hilbert lifting and exact all-momentum weight embeddings yield the complete2e95 M Z136 bound, canonically8e-705. The original two-leg sharp mask with complete actual conversion and lower-band terms gives an anchored error3e54 M Z136/K, canonically12e-746/K. Two inherited matching-table display labels are explicitly corrected without editing their frozen evidence. Independent full complex-dimensional endpoint remainders, frequencies, angular/radial envelopes, noncommuting generator/contact/tangent matrices and norm/omission controls check the result. The bound loses thirteen time and six spatial derivatives. No full lapse/shift/clock/scalar response, reduced inverse, finite-amplitude spatial C2 theorem, finite-coupling quantum background/stability, remaining parent loops/cutoff, finite-gravity IR/Regge or original V/G/B and P8 closure is supplied.",
        "not_established": [
            "Full lapse/shift/clock/scalar Gaussian response and genuinely reduced mixed inverse on appropriate spaces",
            "Finite-amplitude inhomogeneous C2 response, finite-coupling interacting background/stability or physical EFT cutoff",
            "An unanchored homogeneous finite-K rate, strong stress-operator differentiability, global-time bound or same-space contraction",
            "Permission to change the original state, mask, finite prescription, evanescent terms or source germ, or count old/local terms twice",
            "Remaining parent sectors/measure/loops, vacuum matching/cuts/contour/truncation, finite-gravity IR/Regge or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact full dimensional WKB recurrence and explicit complex coefficient bounds, continuous angular/low/near/far domination, common covariance-tangent anchoring, complex tensor/Fourier norm lift and exact regulator subtraction/error. Independent full modes at complex dimension, two-radius endpoint remainders, angular/radial quadrature, literal noncommuting generator/contact/flow tangent identities and omission controls supplement the written proofs. Native/direct/ordinary/CLI retain original SymPy; full regression alone uses the audited exact-GCD adapter. Continuous Fourier, dimensional-limit and functional-analytic arguments are not FORMALIZED. The result remains the tracefree spatial Gaussian sector, not original P8 closure.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete matched tracefree spatial current report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.213.COMPLETE_MATCHED_TRACEFREE_SPATIAL_GAUSSIAN_CURRENT_AND_ANCHORED_UNIFORM_REGULATOR_LIMIT replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
