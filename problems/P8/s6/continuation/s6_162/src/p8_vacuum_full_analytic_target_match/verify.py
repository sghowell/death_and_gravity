"""Read-only full finite-kappa analytic-target classical matching certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_two_loop_elastic_cut import verify as parent

from . import audit, bounds, calibration, holomorphic, jets, transport

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-full-analytic-target-match.json"
PARENT_SHA = "2ef5a0230f712fc5f954ed970e67cef74b86b4216be659075a07b414357d5f13"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_full_analytic_target_match/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen full analytic target parent report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_161_fully_rebuilt": PARENT_SHA,
        "same_finite_kappa_S6_109_rational_step_target_and_S6_111_common_class": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
        "classical_flat_match_not_quantum_target_gravity_or_rolling_parent": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A full analytic target matching proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.162.FULL_FINITE_KAPPA_ANALYTIC_TARGET_CLASSICAL_VACUUM_MATCH",
        "date": "2026-09-10",
        "status": "FULL_FINITE_KAPPA_ANALYTIC_TARGET_CLASSICAL_ACTION_MATCH_ON_RESTRICTED_FLAT_SCHWARTZ_CLASS; NOT_QUANTUM_TARGET_PHYSICAL_TRUNCATION_GLOBAL_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/target.md",
            "notes/jets.md",
            "notes/holomorphic.md",
            "notes/remainder.md",
            "notes/matching.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "full_analytic_target_jets_and_common_holomorphic_circle": serialize(
            {"jets": payload(jets.data()), "holomorphic": payload(holomorphic.data())}
        ),
        "all_higher_field_and_integrated_target_error": serialize(
            {"bounds": payload(bounds.data()), "actual": payload(calibration.data())}
        ),
        "full_classical_target_match_on_same_common_class": serialize(
            payload(transport.data())
        ),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete finite-kappa S6.109 rank-regular rational-step scalar target, not just its quartic germ, is classically matched to the S6.110 polynomial stationary action after the S6.111 field map on the same flat Schwartz/Fourier unit class. Exact sixth/eighth density germs and a uniform complex field-amplitude Cauchy circle bound every remaining even field degree. The target's higher-field action correction is below 1e-1591 times ||Psi||_2^2; adding it to the established field-map/full-resolvent error preserves the strict combined bound below 1e-800 times ||Psi||_2^2. This is not quantum target matching, a rolling-state/cutoff/gravity dictionary or original P8 closure.",
        "not_established": [
            "Quantum matching of the complete analytic target or a physical omitted-loop remainder",
            "A momentum cutoff inferred from the complex field-amplitude circle",
            "Metric variations, stress tensors or finite-gravity G from a fixed flat density estimate",
            "A global inverse map, common rolling-parent B or original P8 closure",
        ],
        "verification_boundary": "Exact canonical weighted germs, coefficient majorants, analytic norm inequalities for the full rational-step target, parity/Cauchy tail and a common-class Parseval bound composed with the unchanged full heavy resolvent. Independent tests derive the actual full-function Taylor coefficients at finite switch orders, directly evaluate the calibrated target at high precision and verify Lorentz invariant and integrated-error controls. Source-pinned written analytic arguments and exact replay are not formalization or independent peer review. Native, direct science, ordinary and CLI retain unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The full analytic target matching report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.162.FULL_FINITE_KAPPA_ANALYTIC_TARGET_CLASSICAL_VACUUM_MATCH replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
