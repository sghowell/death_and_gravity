"""Read-only restricted, tree-propagated CD tensor-noise response bound."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_proca_gaussian import verify as gaussian
from p8_vacuum_affine_tensor_noise_response import verify as parent

from . import audit, control, hamiltonian, response, vertices

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-shear-response.json"
PARENT_SHA = "e815d07a9f855da26c1e26031c55ea040bb521b2bff9c9bed06e1cba7b61f1e5"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_shear_response/*.py"))
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
        "S6_187_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("An actual homogeneous shear-response gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.188.ACTUAL_CD_HOMOGENEOUS_SHEAR_PROCA_HAMILTONIAN_FULL_METRIC_CONTACT_AND_MODEWISE_FINITE_AMPLITUDE_COVARIANCE_REMAINDER",
        "date": "2026-09-11",
        "status": "ACTUAL_SHEAR_HAMILTONIAN_METRIC_CONTACT_AND_MODEWISE_FINITE_AMPLITUDE_COVARIANCE_REMAINDER; NOT_CONTINUUM_RENORMALIZED_RESPONSE_SELF_ENERGY_INVERSE_QUANTUM_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/hamiltonian.md",
            "notes/vertices.md",
            "notes/response.md",
            "notes/remainder.md",
            "notes/state.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_Hamiltonian_and_metric_vertices": serialize(
            {
                "Hamiltonian": payload(hamiltonian.data()),
                "vertices": payload(vertices.data()),
            }
        ),
        "actual_state_response_and_finite_amplitude": serialize(
            {"response": payload(response.data()), "control": payload(control.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged parent on a prescribed homogeneous unimodular shear history has an exact positive three-mode Proca Hamiltonian, with its temporal constraint retained. All metric exponential vertices are energy-normalized and the covariance tangent agrees with the retarded quadratic commutator, including the generally nonzero second-metric contact. The actual background energy propagator is below4 uniformly in momentum; its interaction picture yields an explicit finite-amplitude covariance Taylor remainder for every finite mode in the same all-order Cauchy state. A stated small-shear mode benchmark has relative interaction-picture remainder below1e-6. The general majorant grows in momentum and is not a justified continuum renormalized stress-response bound. Full self-energy feedback, physical cutoff and original V/G/B remain OPEN.",
        "not_established": [
            "A uniform integrable UV-subtracted tensor response, complete covariant finite metric contacts or regulator-independent continuum response remainder",
            "A full quantum metric state, tensor self-energy inverse, scalar/mixed/spatial/nonlinear feedback or controlled interacting background",
            "A physical cutoff from the illustrative mode benchmark, or a dynamic gap from an instantaneous Hamiltonian identity",
            "Full parent measure/loops and matching, finite-gravity IR/Regge or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact constrained Hamiltonian, exponential vertices, physical stress sign, covariance/commutator/contact algebra and continuous energy/Dyson remainder bounds. Independent anisotropic Legendre, rotated polarization, matrix exponential, actual Cauchy-state, finite-amplitude and tangent/contact fixtures supplement the proofs. Native, direct science, ordinary and CLI use original SymPy; only full regression uses the audited exact GCD adapter. Analytic and operator arguments are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The actual homogeneous shear-response report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.188.ACTUAL_CD_HOMOGENEOUS_SHEAR_PROCA_HAMILTONIAN_FULL_METRIC_CONTACT_AND_MODEWISE_FINITE_AMPLITUDE_COVARIANCE_REMAINDER replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
