"""Read-only theorem-applicability and uniform one-particle state certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_flat_dirac_production import verify as parent

from . import audit, scattering, symbols

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-flat-dirac-hadamard.json"
PARENT_SHA = "b9b4c42ca2a054b86797a78e4e6b3e53658288cee87e9cf2b65579015cccd6d0"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_flat_dirac_hadamard/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen free Dirac state parent report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_165_fully_rebuilt": PARENT_SHA,
        "same_specified_free_flat_operator_and_in_out_projectors": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A flat Dirac Hadamard application proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.166.HADAMARD_FREE_FLAT_DIRAC_STATE_AND_UNIFORM_OPERATOR_NORM_IN_OUT_LIMIT",
        "date": "2026-09-10",
        "status": "IDENTIFIED_EXTERNAL_HADAMARD_THEOREM_APPLIES_TO_THE_SPECIFIED_FREE_FLAT_SAT8_DIRAC_IN_OUT_STATES_WITH_UNIFORM_FIRST_QUANTIZED_NORM_LIMIT; NOT_NUMERICAL_STRESS_INTERACTING_CURVED_STATE_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/tails.md",
            "notes/operator.md",
            "notes/theorem.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "all_order_mass_symbol_bounds": serialize(payload(symbols.data())),
        "first_quantized_operator_norm_in_out_limit": serialize(
            payload(scattering.data())
        ),
        "external_hadamard_theorem_application": audit.application(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "All-order Cauchy and compact derivative arguments establish S^-8 mass tails. A uniform bounded-perturbation estimate gives norm-convergent one-particle Moller operators and asymptotic projectors, identifying the same in/out state as the external theorem application. The actual unit-time operator error is below 1e-604 and projector error below 2e-604. The stated free state is Hadamard by the identified external result, not by the finite-energy estimate. No numerical local stress, physical loop error or interacting curved common-parent state is bounded.",
        "not_established": [
            "A numerical bound on the smooth Hadamard remainder or its stress-relevant derivatives",
            "A specified finite local stress-renormalization prescription and quantitative backreaction",
            "A global infinite-volume Fock implementer or interacting gauge asymptotic particle state",
            "The curved bounce's common-parent state, cutoff, physical loop errors, V/G/B or original P8 closure",
        ],
        "verification_boundary": "Exact finite recurrence and Clifford-dictionary checks accompany written all-order analytic and operator-norm arguments. Independent high-order real/complex derivative checks and finite-interval unitary evolution diagnostics test the bounds. The external Hadamard theorem is cited and its hypotheses are matched, not re-proved or formalized. Numerical diagnostics are not validated integration. Native, direct science, ordinary and CLI retain unmodified SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The flat Dirac Hadamard application report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.166.HADAMARD_FREE_FLAT_DIRAC_STATE_AND_UNIFORM_OPERATOR_NORM_IN_OUT_LIMIT replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
