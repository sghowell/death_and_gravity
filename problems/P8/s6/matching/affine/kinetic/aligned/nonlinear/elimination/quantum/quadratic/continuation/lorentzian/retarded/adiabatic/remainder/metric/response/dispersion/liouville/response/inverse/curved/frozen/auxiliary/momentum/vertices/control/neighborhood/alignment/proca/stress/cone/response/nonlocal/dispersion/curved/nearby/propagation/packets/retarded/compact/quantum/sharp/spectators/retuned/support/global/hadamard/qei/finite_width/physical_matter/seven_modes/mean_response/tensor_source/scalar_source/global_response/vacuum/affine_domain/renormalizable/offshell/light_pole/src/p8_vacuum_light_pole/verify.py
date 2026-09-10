"""Read-only actual one-loop polynomial-vacuum light-pole certificate."""

import argparse
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as certificate
from p8_offshell_vacuum import verify as parent

from . import audit, calibration, kernel, normalization

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-one-loop-light-pole.json"
PARENT_SHA = "c471f24ce86d7c8ef9eed0453ecf145fc89efc3b42466c5851a8087cafa66941"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_light_pole/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {key: value for key, value in d.items() if key != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen off-shell vacuum report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_111_fully_rebuilt": PARENT_SHA,
        "same_S6_110_actual_polynomial_model_and_full_Hessian_rebuilt": True,
        "classical_off_shell_field_map_not_used_to_transfer_quantum_measure": True,
        "no_old_cosmological_state_or_counterterm_transfer": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = {key: bool(value) for key, value in audit.gates().items()}
    if not all(value is True for value in gates.values()):
        raise ValueError("An actual one-loop light-pole gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.112.POLYNOMIAL_VACUUM_ONE_LOOP_LIGHT_POLE",
        "date": "2026-09-09",
        "status": "COMPLETE_ACTUAL_ONE_LOOP_LIGHT_TWO_POINT_HESSIAN_NORMALIZATION_FIXED_ON_SHELL_POLE_RESIDUE_ANALYTIC_DISC_AND_POSITIVE_CURVATURE_SCHEME_CONVERSION; NOT_FULL_FOUR_POINT_HIGHER_LOOP_FINITE_GRAVITY_COMMON_BOUNCE_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/normalization.md",
            "notes/analytic.md",
            "notes/subtraction.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": sum(
            value.rows * value.cols if isinstance(value, sp.MatrixBase) else 1
            for value in rows.values()
        ),
        "proof_checks": gates,
        "actual_reduced_Hessian_normalization": serialize(
            payload(normalization.data())
        ),
        "actual_one_loop_light_kernel": serialize(payload(kernel.data())),
        "exact_disc_calibrations": serialize(
            {
                "unit": calibration.point(),
                "half": calibration.point(Fraction(1, 2)),
                "radius_two": calibration.point(2),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The actual full reduced Hessian yields the mixed light-heavy one-loop self-energy with the correct half-trace factor. One fixed on-shell subtraction gives mass one and residue one for the truncated inverse, and a strict complex unit-disc bound excludes additional zeros. The finite local kinetic subtraction is below 3 times 10^-208. Its consistent potential-scheme conversion leaves curvature 1-p with 0<p<10^-405 and preserves the fixed positive quartic and one-loop Phi6 remainder. Complete four-point, higher-loop, finite-gravity and common-bounce matching remain open.",
        "not_established": [
            "All-order propagator poles, residues or a cut-free extension of the one-loop analytic disc",
            "Complete renormalized four-point amplitude, absorptive cuts, contour boundedness or higher-loop error control",
            "Exact stable heavy asymptotic particle or equality with the earlier exponential kinetic parent",
            "Quantum field-map, source, state or counterterm transfer to the affine rolling bounce",
            "Finite-gravity IR/Regge remainder or a controlled common propagating bounce parent",
            "Original P8(b) classification or original P8 closure",
        ],
        "verification_boundary": "Native full reduced-Hessian normalization, literal four-dimensional Feynman shift, anchored analytic kernel identities, convergent radial/parameter integrals and exact rational remainder bounds. Continuum continuation, logarithm-disc and constant-field perturbative subtraction arguments are written/source-pinned, not proof-assistant formalized or peer reviewed. Every ancestor, own source and report field is read-only rebuilt with unmodified scientific SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The polynomial-vacuum light-pole report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.112.POLYNOMIAL_VACUUM_ONE_LOOP_LIGHT_POLE replay passed; full V/G/B and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
