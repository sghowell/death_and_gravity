"""Read-only complete original UV-symbol conversion and all-P remainder."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_proca_gaussian import verify as gaussian
from p8_vacuum_affine_uniform_uv_remainder import verify as parent

from . import angular, audit, centered, flat, geometry

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-subleading-band-conversion.json"
)
PARENT_SHA = "863eeebc61071e48148546d33e778e5af2d7914f51b87688856aae41dcd52aba"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_subleading_band_conversion/*.py"))
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
        "S6_208_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A complete original conversion gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.209.COMPLETE_ORIGINAL_UV_SYMBOL_SHAPE_CONVERSION_ACTUAL_CENTERED_CANCELLATION_AND_UNIFORM_TAIL",
        "date": "2026-09-11",
        "status": "COMPLETE_ORIGINAL_UV_SYMBOL_SHAPE_CONVERSION_ACTUAL_CENTERED_CANCELLATION_AND_UNIFORM_TAIL; NOT_FULL_CURVED_A1_EVALUATION_CONTACT_COVARIANT_MATCHING_RESPONSE_INVERSE_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/shell.md",
            "notes/angular.md",
            "notes/symmetry.md",
            "notes/benchmark.md",
            "notes/remainder.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_original_shell_angular_shapes_and_full_massive_flat_benchmark": serialize(
            {
                "geometry": payload(geometry.data()),
                "angular": payload(angular.data()),
                "flat": payload(flat.data()),
            }
        ),
        "actual_centered_cancellation_and_uniform_original_conversion_tail": serialize(
            payload(centered.data())
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The exact original two-created-mode shell retains ten universal shapes, twenty j/d/h branches and thirty-five source time-jet shapes through finite cutoff order. The positive grazing strip contributes its cubic value and quartic derivative. The actual full unit-W8 pair product obeys exchange and Schwarz symmetries with independent detector/source times and distinct sector phases. In centered momenta these imply even total ultraviolet grades, not deletion of odd endpoint labels. The complete original UV-symbol conversion consequently has cubic and linear cutoff terms but no quadratic or finite term. All individual terms are first retained; a nonzero leading-slot finite artifact is not substituted for the full sum. An explicit angular Cauchy/strip bound, degree-eight Chebyshev control and complementary large-transfer power/log estimates give a uniform conversion remainder1e40||D||L2 X46[Gamma]/K for every external momentum and K>=2m. Full all-angle massive flat fields independently check all45 ultraviolet coefficients and explicit nonzero cubic/linear channel terms. Actual curved full-field/projector, exchange/parity and two-resolution ultraviolet tests retain additional allowed-slot cancellations and nonzero finite-momentum odd endpoints. The full curved linear coefficient remains its actual centered f2c functional and still needs evaluation. The complete one-leg contact, one-ball UV finite/divergent coefficients and original fixed covariant matching are not replaced or completed. No full matched response, mixed inverse, finite-coupling background/stability, physical cutoff, remaining parent loops, finite-gravity IR/Regge or original V/G/B closure follows. Original P8 remains open.",
        "not_established": [
            "Explicit full curved linear cutoff-conversion coefficient, complete one-ball UV/contact finite and divergent coefficients and original fixed covariant matching",
            "Permission to replace the original regulator, alter the full W8 comparison or actual state, change finite counterterms or double-count earlier finite pieces",
            "A full matched response, reduced lapse/shift or derivative-compatible mixed inverse, finite-coupling background or stability",
            "Remaining parent loops/cutoff/scattering/IR/Regge or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact shell/grazing moments, angular filtration identities, all45 flat coefficients, full generic exchange/Schwarz identities, centered cancellation, finite Chebyshev derivative bounds and explicit all-P power/log tail constants. Independent original shell/cap quadratures, all-angle complete massive flat fields, actual full curved modes and two-resolution coefficients, and every original UV slot in small/large transfer regimes supplement written continuous arguments. Native, direct science, ordinary and CLI retain original SymPy; only full regression uses the audited exact GCD adapter. Continuous analyticity, angular/Fourier norm and regulator removal are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The original complete conversion report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.209.COMPLETE_ORIGINAL_UV_SYMBOL_SHAPE_CONVERSION_ACTUAL_CENTERED_CANCELLATION_AND_UNIFORM_TAIL replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
