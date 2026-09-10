"""Read-only polynomial vacuum and positive one-loop potential certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as certificate
from p8_affine_vacuum_domain import verify as parent

from . import audit, calibration, gaussian, model, potential, powercount

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-quantum-preparation.json"
PARENT_SHA = "a589a7b8a48bb0bfda6b5a4aa36e9262a3006ad3a1faaba32c311e7848992eaf"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_polynomial_vacuum/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {key: value for key, value in d.items() if key != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen regular affine vacuum-domain report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_109_fully_rebuilt": PARENT_SHA,
        "S6_108_actual_contact_and_exact_full_heavy_exchange_rebuilt": True,
        "new_polynomial_vacuum_not_identified_with_old_affine_or_kinetic_parent": True,
        "frozen_action_state_and_counterterm_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = {key: bool(value) for key, value in audit.gates().items()}
    if not all(value is True for value in gates.values()):
        raise ValueError("A polynomial vacuum or quantum preparation gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.110.POLYNOMIAL_VACUUM_QUANTUM_PREPARATION",
        "date": "2026-09-09",
        "status": "COMPLETE_NEW_STABLE_CANONICAL_POLYNOMIAL_VACUUM_EXACT_FULL_ON_SHELL_TREE_FOUR_POINT_MATCH_FINITE_REGULATOR_GAUSSIAN_REDUCTION_COUNTERTERM_POWERCOUNT_AND_POSITIVE_BOUNDED_ONE_LOOP_CONSTANT_POTENTIAL; NOT_FULL_LOOP_AMPLITUDE_V_G_B_COMMON_BOUNCE_PARENT_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/tree.md",
            "notes/gaussian.md",
            "notes/powercount.md",
            "notes/potential.md",
            "notes/scope.md",
            "notes/literature.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": sum(
            v.rows * v.cols if isinstance(v, sp.MatrixBase) else 1
            for v in rows.values()
        ),
        "proof_checks": gates,
        "polynomial_model": serialize(payload(model.data())),
        "exact_finite_regulator_integration": serialize(payload(gaussian.data())),
        "power_counting_counterterm_basis": serialize(payload(powercount.data())),
        "positive_one_loop_potential": serialize(payload(potential.data())),
        "exact_calibrations": serialize(
            {
                "actual_parameters": calibration.point(),
                "independent_rational_parameters": calibration.point(2, 1),
                "zero_momentum": calibration.kernel_point(),
                "zero_field": calibration.kernel_point(1, 0),
                "large_Euclidean_momentum": calibration.kernel_point(
                    10**197, -sp.Rational(7, 3)
                ),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "A distinct globally stable canonical polynomial massive vacuum realizes the full earlier on-shell tree exchange. Exact Gaussian heavy integration and the complete Schur Hessian retain all mixed quantum contributions in a nonlocal light functional. Superficial graph counting preserves the Gaussian-heavy counterterm structure. A fixed-subtraction constant-background one-loop remainder is positive and bounded by an explicit rational Phi-sixth-power ceiling. None of these supplies the full renormalized scattering/dispersion error, finite-gravity contour or common propagating bounce parent.",
        "not_established": [
            "Quantum on-shell light pole mass/residue, complete renormalized b2 correction, absorptive/cut and higher-loop/truncation bounds",
            "Holomorphic loop extension of the tree polydisc, an exactly stable heavy asymptotic particle or a bare-pole-squared dispersion integral",
            "All-order perturbation convergence or nonperturbative continuum all-energy UV completion",
            "General inhomogeneous Hessian positivity or uniform relative loop control at arbitrarily large fields",
            "Identity of the polynomial vacuum, derivative-exchange vacuum and regular auxiliary affine bounce parent beyond the stated tree contact",
            "Finite-gravity IR/Regge remainder, controlled common bounce matching, old quantum state/counterterm transfer or original P8(b)/P8 closure",
        ],
        "verification_boundary": "Native symbolic tree amplitude, literal potential derivatives, independent two-site Gaussian action/Hessian and 24-permutation determinant, exact graph-incidence counting, logarithmic integral identities and rational calibrations. Continuous positivity, regulator, perturbative power-counting and integral arguments are written/source-pinned, not proof-assistant formalized or peer reviewed. Complete read-only ancestor, own-source and report-field replay uses unmodified scientific SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The polynomial vacuum quantum-preparation report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.110.POLYNOMIAL_VACUUM_QUANTUM_PREPARATION replay passed; full V/G/B and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
