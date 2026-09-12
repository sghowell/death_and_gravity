"""Read-only full spatial flat covariant/time-subtraction conversion."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_flat_tensor_cut import verify as parent
from p8_vacuum_affine_proca_gaussian import verify as gaussian

from . import audit, bounds, conversion, flat, tensor

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-flat-spatial-conversion.json"
PARENT_SHA = "e7a37017276d0b5474d81a295ab3f2574f69b883ffee34b032a0dadf9566ceec"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_flat_spatial_conversion/*.py"))
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
        "S6_199_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A full spatial flat subtraction conversion gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.200.FULL_SPATIAL_FLAT_PROCA_COVARIANT_AND_RETARDED_TIME_SUBTRACTION_CONVERSION_WITH_NONLOCAL_COEFFICIENT_BOUNDS",
        "date": "2026-09-11",
        "status": "FULL_SPATIAL_FLAT_PROCA_COVARIANT_AND_RETARDED_TIME_SUBTRACTION_CONVERSION_WITH_NONLOCAL_COEFFICIENT_BOUNDS; NOT_CURVED_SPATIAL_MATCHING_PHYSICAL_LOCAL_POLYNOMIAL_FULL_INVERSE_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/tensor.md",
            "notes/conversion.md",
            "notes/flat.md",
            "notes/bounds.md",
            "notes/limit.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_tensor_and_exact_subtraction_conversion": serialize(
            {"tensor": payload(tensor.data()), "conversion": payload(conversion.data())}
        ),
        "retarded_origin_and_uniform_spatial_bounds": serialize(
            {"flat": payload(flat.data()), "bounds": payload(bounds.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The full mass1000 canonical Proca flat spatial cut admits polynomial covariant tensor numerators with no lightcone projector pole. Exact cubic division expresses the specified covariant nonlocal representative as the actual flat sixth-time-derivative subtracted spectral bulk plus three convergent equal-time tensor coefficients. Their nonzero-transfer dependence is genuinely spatially nonlocal, with all-momentum coefficient bounds and explicit inverse-Lambda spectral tails. The standard six-spatial/four-time source norm gives conversion bound1e10 and tail1e17/Lambda; both external-metric canonical factors yield4e-790 and4e-783/Lambda. Combining this with the absolute spectral sine bulk gives the complete specified flat nonlocal weak bound2e10 in the stated six-spatial/six-time source norm, with tail1e14/sqrt(Lambda) and canonical displays8e-790 and4e-786/sqrt(Lambda). Complete nonzero-transfer curvature Hessians, matrix spectral integrals, frequency derivatives and prepared sine-kernel tests are independent checks. This representative conversion does not select the physical local polynomial, identify different finite regulators, remove curved odd endpoints, solve the actual CD covariant contact/endpoint matching, provide a mixed same-space inverse, finite-coupling background, stability, physical cutoff or finite-gravity Regge data. Original V/G/B and P8 remain open.",
        "not_established": [
            "The physical finite local tensor/contact polynomial or full curved CD reference contact/endpoint matching",
            "Identification of the original one-mode contact and two-mode memory finite regulators",
            "A full mixed response inverse, finite-coupling interacting background, stability or physical cutoff matching",
            "Full parent scattering cuts, contour/truncation, finite-gravity IR/Regge or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact complete spatial polynomial numerators, generic cubic division, spectral residues, moment constants, canonical factors and retarded identities. Independent covariant four-curvature tensors at nonzero transfer, full matrix improper spectral integrals, centered frequency derivatives, operator norms, infinite tails and literal prepared sine kernels supplement written analytic arguments. Native, direct science, ordinary and CLI use original SymPy; only full regression uses the audited exact GCD adapter. The continuous analytic and operator-norm arguments are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The full spatial flat subtraction conversion report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.200.FULL_SPATIAL_FLAT_PROCA_COVARIANT_AND_RETARDED_TIME_SUBTRACTION_CONVERSION_WITH_NONLOCAL_COEFFICIENT_BOUNDS replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
