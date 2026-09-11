"""Read-only actual matrix adiabatic reference and uniform covariance tail."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_local_tensor_response import verify as parent
from p8_vacuum_affine_proca_gaussian import verify as gaussian

from . import audit, frame, initial, jets, riccati

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-matrix-adiabatic.json"
PARENT_SHA = "567a8a84e981c637b000ef671c350b63c54f4fc357ca60eaf1ce1de9c1293b02"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_matrix_adiabatic/*.py"))
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
        "S6_189_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("An actual matrix adiabatic covariance-tail gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.190.ACTUAL_CD_FULL_MATRIX_ADIABATIC_REFERENCE_AND_UNIFORM_ALL_ORDER_STATE_HIGH_MOMENTUM_COVARIANCE_REMAINDER",
        "date": "2026-09-11",
        "status": "ACTUAL_UNCHANGED_STATE_FULL_MATRIX_ADIABATIC_HIGH_BAND_COVARIANCE_REMAINDER; NOT_RESPONSE_PARAMETER_REMAINDER_COVARIANT_SUBTRACTION_INVERSE_FINITE_COUPLING_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/frame.md",
            "notes/jets.md",
            "notes/reference.md",
            "notes/state.md",
            "notes/error.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "full_matrix_frame_and_uniform_jets": serialize(
            {"frame": payload(frame.data()), "jets": payload(jets.data())}
        ),
        "finite_reference_and_actual_state_tail": serialize(
            {"reference": payload(riccati.data()), "state": payload(initial.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The actual constrained three-mode Hamiltonian admits a momentum-uniform matrix canonical frame and finite tenth-order Riccati reference retaining noncommuting polarization mixing. Explicit square-root Sylvester and finite-jet bounds avoid a condition-number or polarization-gap assumption. Quantitative matching retains the original all-order Cauchy state and all higher cutoff-weighted terms. On the specified homogeneous C12 shear neighborhood, the actual full balanced covariance differs from the reference by less than2e31 nu_minus^-10 for nu_minus>=1e16; the complete infinite energy-weighted covariance tail is below1e-65. This is not yet an amplitude-parameter response remainder, fixed-covariant subtraction/contact match, full feedback inverse, finite-coupling quantum background or physical cutoff. Original V/G/B and P8 remain OPEN.",
        "not_established": [
            "Amplitude-parameter response derivatives or their finite-coupling Taylor remainder",
            "The complete fixed-covariant UV subtraction and remaining finite response contacts",
            "Low-band response, full spatial/mixed/nonlinear feedback, inverse or interacting quantum background",
            "Physical cutoff/threshold, full interacting parent matching, finite-gravity IR/Regge or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact full matrix flow and pure-state covariance, finite real-time jet and ordered Riccati recurrences, actual all-order initial matching and continuous infinite radial estimates. Independent ill-conditioned frame, noncommuting covariance/reference and actual T/L preparation fixtures supplement the written proofs. Native, direct science, ordinary and CLI use original SymPy; only full regression uses the audited exact GCD adapter. Analytic and functional arguments are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The actual matrix adiabatic covariance report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.190.ACTUAL_CD_FULL_MATRIX_ADIABATIC_REFERENCE_AND_UNIFORM_ALL_ORDER_STATE_HIGH_MOMENTUM_COVARIANCE_REMAINDER replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
