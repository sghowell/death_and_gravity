"""Read-only explicit local finite quantum regulator with unchanged ancestry."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_nonlinear_spatial_reduction import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-local-quantum-regulator.json"
PARENT_SHA = "36be009c0d58f8e99a4b703bf5a014dde6a28c07dbb888453a8d04aa0fba7f98"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_local_quantum_regulator/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen full nonlinear spatial reduction report changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_267_entire_spatial_cotangent_reference_report_and_ancestry_rebuilt": PARENT_SHA,
        "original_source_parameters_and_pure_nonzero_preparations_unchanged": True,
        "S261_refutation_and_S265_Gaussian_integrability_boundary_unchanged": True,
        "new_finite_quantum_regulator_not_original_continuum_matching": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "A full local quantum regulator or unchanged-reference gate failed"
        )
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.268.COMPLETE_EXPLICIT_FINITE_LOCAL_SOURCE_PINNED_QUANTUM_REGULATOR_WITH_FULL_CANONICAL_TIME_AND_WEYL_ORDERING_CONTACTS_SAME_FREE_REFERENCE_EXACT_UNITARY_DYNAMICS_AND_DIMENSION_DEPENDENT_STATE_COMPARISON_CRITERIA_WITH_UNCHANGED_PARENT",
        "date": "2026-09-14",
        "status": "EXPLICIT_FINITE_LOCAL_REGULATED_UNITARY_QUANTUM_DYNAMICS_AND_EXACT_ERROR_CRITERIA; NOT_EVALUATED_PHYSICAL_LEAKAGE_ORDERING_MATCHING_REGULATOR_REMOVAL_OR_ORIGINAL_V_G_B_P8",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/ordering.md",
            "notes/localization.md",
            "notes/dynamics.md",
            "notes/reference.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_source_time_connection_full_ordering_and_physical_contacts": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_phase_localization_and_convergent_quantum_dynamics": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "An explicitly named finite IR/UV/phase-localized quantization retains the entire original Hamiltonian, nonlinear auxiliary/spatial reconstruction, canonical time connection, original free flow and same pure nonzero seed. Full-covariance coherent quantization with declared first Weyl calibration produces bounded self-adjoint interactions and actual unitary dynamics preserving all residual translations. The positive finite coherent volume and the separately calibrated readout are distinguished. Exact ordering, dimension-dependent localization, Dyson and state-comparison criteria retain every contact. No actual P8 domain radius, time width, small leakage, ordering size or original interacting mean is evaluated. The regulator changes continuation outside the local chart and does not define original quantum matching or a continuum limit. Original V/G/B/P8 remain OPEN.",
        "not_established": [
            "An evaluated common canonical/whitened radius or physical time interval",
            "Actual small P8 state leakage, interaction norms or ordering error",
            "Unconditional positivity of the first-Weyl-calibrated volume",
            "Equality with an original unregularized singular Hamiltonian or its quantum mean",
            "A homogeneous quantum state or identity of the finite family with the R3 state",
            "A continuum diffeomorphism/BRST regulator, regulator removal or loop estimates",
            "Original Wilsonian matching, finite-gravity gates or nonlinear completeness",
            "Formal verification of the written functional and operator arguments",
        ],
        "verification_boundary": "Exact whole-parent/time bindings, full-covariance heat and physical Hessian contacts, complete dimension-dependent tail identities and noncommuting coherent/Dyson fixtures support independent written proofs. Gaussian quadrature and matrix ODE comparisons are diagnostics, not full P8 dynamics or a finite CCR representation. Native/direct/ordinary/CLI use original SymPy; only the captured full regression uses its audited exact-GCD adapter. Every frozen ancestor, historical erratum and original open gate remains unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete local quantum regulator report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.268 explicit finite local quantum regulator replay passed; quantitative original matching and P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
