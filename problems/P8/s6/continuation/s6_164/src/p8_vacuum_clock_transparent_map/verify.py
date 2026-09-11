"""Read-only polynomial clock-transparent map and common-class certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_protected_yukawa_profile import verify as parent

from . import audit, calibration, counting, gate, norms, transport

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-clock-transparent-map.json"
PARENT_SHA = "56a114490c2334adcdfbfb14e1a2add42e7f7e769e51b024ff91141ec63074cb"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_clock_transparent_map/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen clock-transparent map parent report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_163_fully_rebuilt": PARENT_SHA,
        "same_SAT8_lower_field_reference_boundary_and_mass_profile": True,
        "new_polynomial_source_coordinate_not_new_parent_bounce_solution": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A clock-transparent map proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.164.CLOCK_TRANSPARENT_POLYNOMIAL_MAP_AND_FULL_TARGET_CLASSICAL_MATCH",
        "date": "2026-09-10",
        "status": "CLOCK_TRANSPARENT_POLYNOMIAL_PHYSICAL_SOURCE_MAP_WITH_SAME_NAMED_TWO_LOOP_VACUUM_DATA_AND_CONTROLLED_FULL_TARGET_CLASSICAL_MATCH; NOT_BACKGROUND_SOLUTION_QUANTUM_STATE_PHYSICAL_REMAINDER_GLOBAL_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/gate.md",
            "notes/grading.md",
            "notes/norms.md",
            "notes/matching.md",
            "notes/clock.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "polynomial_gate_and_two_domain_grading": serialize(
            {"gate": payload(gate.data()), "counting": payload(counting.data())}
        ),
        "common_class_full_resolvent_action_comparison": serialize(
            {"norms": payload(norms.data()), "actual": payload(calibration.data())}
        ),
        "named_vacuum_and_clock_argument_dictionary": serialize(
            payload(transport.data())
        ),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The finite polynomial map Fhat=Psi+q8(X)R_old, with X=(partial Psi)^2/kappa, equals the canonical clock identically at physical X=1 and has seven vanishing nonidentity variations there. Its first vacuum map change has field degree19, so the named SAT8 through-two-loop physical-source coefficients remain unchanged. On the original Schwartz/Fourier unit class its finite highest degree33 gives squared-source radius66, still strictly inside the full heavy gap. An explicit action comparison, including all free/local/heavy terms, bounds the new map error below 1e-6790 times ||Psi||_2^2 and preserves the complete analytic-target match below 1e-800 times ||Psi||_2^2. The protected fermion mass-derivative ratios now apply to the same kinematic unit-clock argument. No parent background equation, quantum rolling state, physical truncation, common-parent B or original P8 closure is inferred.",
        "not_established": [
            "The target bounce as a solution of the ordinary canonical polynomial parent",
            "Global nonlinear invertibility or a full rolling quantum state/particle-production bound",
            "Quantum target matching, physical omitted-order errors or an interacting cutoff",
            "Finite-gravity G, common-parent B, global V contours or original P8 closure",
        ],
        "verification_boundary": "Exact interpolation and clock/vacuum jets, finite polynomial support counting, common-class L2/Linfinity inequalities, unexpanded heavy-resolvent norm and rational actual error composition. Independent tests integrate the interpolation polynomial, retain nonzero eighth clock variations, enumerate finite Fourier convolutions and compare full stationary actions and linearized majorants. Source-pinned written analytic arguments and exact replay are not formalization or independent peer review. Native, direct science, ordinary and CLI retain unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The clock-transparent map report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.164.CLOCK_TRANSPARENT_POLYNOMIAL_MAP_AND_FULL_TARGET_CLASSICAL_MATCH replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
