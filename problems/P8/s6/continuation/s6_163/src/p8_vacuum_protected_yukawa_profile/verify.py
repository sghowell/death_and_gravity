"""Read-only GY14-SAT8 protected mass profile and named coefficient transport."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_full_analytic_target_match import verify as parent

from . import audit, calibration, clock, order, profile, transport

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-protected-yukawa-profile.json"
PARENT_SHA = "50d2439dea8d9b3889f4875679c05f27d2d4e9e81ef88d33bc404cb4671e3cfd"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_protected_yukawa_profile/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen protected profile parent report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_162_fully_rebuilt": PARENT_SHA,
        "same_actual_lower_field_GY14_reference_boundary": True,
        "new_SAT8_EFT_not_original_global_linear_Yukawa_model": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A protected mass profile proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.163.GY14_SAT8_PROTECTED_MASS_PROFILE_AND_NAMED_TWO_LOOP_INVARIANCE",
        "date": "2026-09-10",
        "status": "NEW_SAT8_ANALYTIC_EFT_WITH_GLOBAL_REAL_POINTWISE_MASS_FLOOR_AND_IDENTICAL_NAMED_TWO_LOOP_VACUUM_DATA; NOT_ROLLING_STATE_PHYSICAL_REMAINDER_GLOBAL_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/candidate.md",
            "notes/profile.md",
            "notes/counting.md",
            "notes/reference.md",
            "notes/clock.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "new_profile_and_dimensional_loop_grade": serialize(
            {"profile": payload(profile.data()), "order": payload(order.data())}
        ),
        "actual_mass_floor_and_conditional_clock_screen": serialize(
            {"clock": payload(clock.data()), "actual": payload(calibration.data())}
        ),
        "identical_named_coefficients_and_classical_match": serialize(
            payload(transport.data())
        ),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The separately named GY14-SAT8 EFT replaces the two active opposite linear Yukawa arguments by f_R(Phi)=Phi/(1+(Phi/R)^8)^(1/8), R=10^300. Both real pointwise fermion masses remain strictly between 0.99mF and 1.01mF for every real argument. New vertices begin with nine scalar legs; connected loop and forest grading proves that the named physical Phi pole, four-point/b2, vacuum/H references and low-energy two-loop cut retain exactly their previous through-two-loop coefficients. The classical zero-fermion full-target match is unchanged. A conditional direct-argument clock bound is not transferable to the literal cubic-map extrapolation: its flat linear-clock first mass-derivative ratio exceeds 1e97 at zero. Neither this smallness-screen failure nor the mass floor establishes a state, physical remainder, UV completion, common-parent B or original P8 closure.",
        "not_established": [
            "A global Dirac propagator gap, transported rolling state or particle-production bound",
            "Equality of arbitrary higher-point or higher-loop functions with original GY14",
            "Global renormalizability, inherited all-field UV running or a physical interacting cutoff for SAT8",
            "Quantum target matching, physical omitted-order errors, finite-gravity G, common-parent B or original P8 closure",
        ],
        "verification_boundary": "Exact real profile derivatives and vacuum binomial germ, loop/forest and dimensional grading, rational calibrated pointwise mass and conditional derivative bounds, and literal flat-jet cubic-map screen. Independent tests evaluate the full profile, extract higher-field determinant changes, enumerate compatible mixed-field halfedges and compare regulated Gaussian perturbations. Source-pinned written analytic arguments and exact replay are not formalization or independent peer review. Native, direct science, ordinary and CLI retain unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The protected mass profile report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.163.GY14_SAT8_PROTECTED_MASS_PROFILE_AND_NAMED_TWO_LOOP_INVARIANCE replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
