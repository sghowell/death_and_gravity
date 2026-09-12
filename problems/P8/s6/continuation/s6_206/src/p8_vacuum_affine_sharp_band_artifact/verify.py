"""Read-only actual leading sharp two-leg endpoint regulator conversion."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_finite_endpoint_polynomial import verify as parent
from p8_vacuum_affine_proca_gaussian import verify as gaussian

from . import angular, asymptotics, audit, geometry, polarizations

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-sharp-band-artifact.json"
PARENT_SHA = "9e5027a0d817d11ae196410e0863f0e4a93d3d647a8caaeb85a37e1cfb65e71e"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_sharp_band_artifact/*.py"))
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
        "S6_205_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("An actual sharp-band artifact gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.206.ACTUAL_FULL_PROCA_LEADING_SHARP_TWO_LEG_ENDPOINT_REGULATOR_CONVERSION",
        "date": "2026-09-11",
        "status": "ACTUAL_FULL_PROCA_LEADING_SHARP_TWO_LEG_ENDPOINT_REGULATOR_CONVERSION; NOT_FULL_REGULATOR_COVARIANT_MATCHING_PHYSICAL_DIVERGENCE_INVERSE_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/polarizations.md",
            "notes/angular.md",
            "notes/geometry.md",
            "notes/asymptotics.md",
            "notes/conversion.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "actual_full_polarization_symbol_and_angular_coefficient": serialize(
            {
                "polarizations": payload(polarizations.data()),
                "angular": payload(angular.data()),
            }
        ),
        "original_sharp_geometry_and_actual_asymptotic_conversion": serialize(
            {
                "geometry": payload(geometry.data()),
                "asymptotics": payload(asymptotics.data()),
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The actual massive Proca unit-W8 first endpoint, summed over all nine physical polarization pairs with the current sign and Fourier measure retained, has leading density r*[4tr(DG)-4n.(DG+GD)n+3(n.D.n)(n.G.n)]/(8a). Its original two-created-mode band differs from an explicitly named one-k-ball comparator by a real leading lost-shell contribution. The first-five-endpoint conversion is -K^3|P|[18tr(DG)-12(D phat).(G phat)-(phat.D.phat)(phat.G.phat)]/(512 pi^2 a)+o(K^3) at fixed nonzero P on the compact CD slab. The bracket eigenvalues on all five normalized tracefree channels are18,18,12,12,28/3. Exact weighted angular moments, grazing geometry and uniform actual W8 asymptotics justify the coefficient; the other four endpoint differences are lower shell order and the original one-leg contact is unchanged by this comparison. The longitudinal contribution is retained: this is not Maxwell theory. The nonzero |P| cusp is a regulator-conversion artifact, not a finite local Hessian. The original regulator is not replaced, no subtraction is added here and the comparator is not claimed covariantly renormalized. This result does not establish all subleading artifacts, finite/divergent quantum matching, a divergence of the renormalized physical current, a model exclusion, full response/inverse/background/stability, cutoff or original V/G/B closure. The canonical prefactor does not make K^3 uniformly small as the regulator is removed. Original P8 remains open.",
        "not_established": [
            "Complete subleading sharp-band artifacts and finite/divergent quantum contact/candidate-cell coefficients in the fixed covariant prescription",
            "A renormalized one-ball comparator, permission to replace the original regulator, a physical current divergence or model exclusion",
            "A uniform all-external-momentum norm, full matched response, derivative-compatible inverse, finite-coupling background or stability",
            "Full parent loops/cutoff/scattering/IR/Regge or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact all-nine-pair leading Proca symbol, complete hemisphere moments, all five spatial channels, original shell geometry, full W8 scaled limits, current sign and canonical factors. Independent actual full W8 modes at three momentum scales, multiple compact times, directions and noncommuting tensors, original lost-shell quadratures, rotations, longitudinal and nonpolynomial controls and prewarmed exact-order guards supplement the written uniform fixed-P asymptotic and dominated shell arguments. Native, direct science, ordinary and CLI use original SymPy; only full regression uses the audited exact GCD adapter. These continuous asymptotic arguments are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The actual sharp-band artifact report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.206.ACTUAL_FULL_PROCA_LEADING_SHARP_TWO_LEG_ENDPOINT_REGULATOR_CONVERSION replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
