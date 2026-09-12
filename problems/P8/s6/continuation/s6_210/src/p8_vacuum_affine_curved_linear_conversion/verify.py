"""Read-only actual complete curved linear original-cutoff coefficient."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_proca_gaussian import verify as gaussian
from p8_vacuum_affine_subleading_band_conversion import verify as parent

from . import angular, audit, bounds, density, jets

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-curved-linear-conversion.json"
)
PARENT_SHA = "0ffbe7d77d194bbae28237b7fe54cda45b6ffe0630205628e3992d5f345e9f0b"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_curved_linear_conversion/*.py"))
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
        "S6_209_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("An actual curved linear coefficient gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.210.ACTUAL_FULL_CURVED_LINEAR_CUTOFF_CONVERSION_COEFFICIENT_ALL_TENSOR_CHANNELS_AND_SOURCE_TIME_JETS",
        "date": "2026-09-11",
        "status": "ACTUAL_FULL_CURVED_LINEAR_CUTOFF_CONVERSION_COEFFICIENT_ALL_TENSOR_CHANNELS_AND_SOURCE_TIME_JETS; NOT_FULL_ONE_BALL_UV_CONTACT_COVARIANT_MATCHING_RESPONSE_INVERSE_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/jets.md",
            "notes/density.md",
            "notes/angular.md",
            "notes/operator.md",
            "notes/bounds.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "actual_complete_WKB_and_curved_second_grade_source_densities": serialize(
            {"jets": payload(jets.data()), "density": payload(density.data())}
        ),
        "full_tensor_curved_conversion_operator_and_coefficient_norms": serialize(
            {"angular": payload(angular.data()), "bounds": payload(bounds.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The actual complete curved linear original-cutoff coefficient is now evaluated, not replaced by the massive flat benchmark. Exact frozen WKB jets retain P1_T(z1)=0 and P1_L(z1)=-(H_prime+2H^2)/2. Full transverse, mixed and longitudinal constrained amplitudes give the actual second-grade curvature correction; the imaginary mixed amplitudes are not discarded. The j1/d1 slot evaluates to zero, while its higher nonzero slot and finite odd endpoints remain. Independent source/detector times give the full j2/d0 first and second source-time derivatives. Seven general-tensor azimuthal projector identities reconstruct all five tracefree channels. The resulting A1 includes mass, spatial, curvature and a*(Gamma_second+H Gamma_prime) terms in four explicit invariant forms. The original A3, zero complete quadratic/finite conversion and all-P remainder1e40 X46/K remain unchanged. Proper-time divergence/Green identities and explicit coefficient bounds1/100 for A3 and1e4 for A1 are established in stated norms, with canonical displays4e-802 and4e-796. These coefficient bounds are not uniform after multiplication by K3 and K. Independent full W8 curved/mass-scaled-flat fields, multiple clocks/directions/noncommuting tensors, source time jets, complete angular quadrature and two Cauchy resolutions check the formula. The full original one-ball UV/contact finite and divergent coefficients and fixed covariant matching still require completion. No full matched response, mixed inverse, finite-coupling background/stability, physical cutoff, remaining parent loops, finite-gravity IR/Regge or original V/G/B closure follows. Original P8 remains open.",
        "not_established": [
            "Complete original one-ball UV/contact finite and divergent coefficients and fixed covariant matching",
            "Permission to replace the original regulator or reference/state, change finite counterterms, or add old finite pieces or the local target twice",
            "A full matched response, reduced mixed inverse, finite-coupling background/stability, or cutoff-uniform smallness from coefficient normalization",
            "Remaining parent loops/cutoff/scattering/IR/Regge or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact actual WKB low jets, full two-time physical amplitudes, complete source-time iterations, seven general-tensor azimuthal reconstructions, all curved invariant coefficients, proper-time Green identities and explicit norm/canonical constants. Independent full-mode coefficient extraction at separated transfer scales, all source time jets, actual angular integrals at two inverse-radius resolutions and exact finite Gauss moments supplement the written continuous estimates. Native, direct science, ordinary and CLI retain original SymPy; only full regression uses the audited exact GCD adapter. Continuous Fourier-norm and regulated-current arguments are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The actual curved linear coefficient report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.210.ACTUAL_FULL_CURVED_LINEAR_CUTOFF_CONVERSION_COEFFICIENT_ALL_TENSOR_CHANNELS_AND_SOURCE_TIME_JETS replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
