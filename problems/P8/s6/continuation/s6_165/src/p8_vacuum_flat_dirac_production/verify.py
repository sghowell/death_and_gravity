"""Read-only free flat-clock Dirac production and transition-tail certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_clock_transparent_map import verify as parent

from . import audit, calibration, energy, mode, profile, transition

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-flat-dirac-production.json"
PARENT_SHA = "2684caa2dd2764c3fcc525c31b571188b323ca8c2cc81da0ce8eb9b2e69bd1b0"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_flat_dirac_production/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen free Dirac production parent report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_164_fully_rebuilt": PARENT_SHA,
        "same_protected_mass_profile_and_aligned_clock_argument": True,
        "specified_quadratic_flat_subsystem_not_full_interacting_state": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A free Dirac production proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.165.FREE_FLAT_CLOCK_DIRAC_OUT_ENERGY_AND_COMPLETE_TRANSITION_REMAINDER",
        "date": "2026-09-10",
        "status": "SPECIFIED_FREE_FLAT_CLOCK_DIRAC_IN_OUT_STATE_WITH_COMPLETE_TRANSITION_SERIES_REMAINDER_AND_FINITE_OUT_PARTICLE_ENERGY; NOT_INTERACTING_OR_CURVED_STATE_LOCAL_STRESS_PHYSICAL_LOOP_REMAINDER_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/state.md",
            "notes/profile.md",
            "notes/transition.md",
            "notes/energy.md",
            "notes/validation.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "free_mode_state_and_switching_profile": serialize(
            {"mode": payload(mode.data()), "profile": payload(profile.data())}
        ),
        "complete_transition_remainder_and_out_energy": serialize(
            {"transition": payload(transition.data()), "energy": payload(energy.data())}
        ),
        "actual_free_clock_energy_and_reference_scale_comparison": serialize(
            payload(calibration.data())
        ),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "For the specified quadratic flat Dirac subsystem with masses mF +/- Delta s(t/tau), the protected SAT8 profile has positive in/out masses and integrable tails, giving modewise in/out evolution. Two integrations by parts bound the first transition integral uniformly in momentum. The complete odd Dyson tail in the transition coupling is retained, not discarded. Both terms decay as p^-3; the exact spin/flavor/particle count and radial integrals give free out-particle energy below 1e784 and a ratio below 1e-16 to the named kappa reference scale. The uniform transition amplitude is below 1e-9. These are exact bounds for the stated quadratic evolution, not all Feynman loops, an interacting gauge asymptotic state, a local renormalized stress bound, a curved bounce state or original P8 closure.",
        "not_established": [
            "Gauge-invariant asymptotic colored particles or the full interacting GY14 state",
            "A curved FLRW production estimate, transient renormalized local stress/coherence or backreaction",
            "A Hadamard remainder or a single global infinite-volume Fock unitary from finite occupation energy",
            "Physical omitted Feynman-loop errors, quantum target matching, V/G/B or original P8 closure",
        ],
        "verification_boundary": "Exact two-level rotation and mode norm, full-profile total-variation estimates, boundary-controlled integrations by parts, convergent odd transition-series tail and complete momentum energy integration with rational calibrated parameters. Independent tests solve finite-interval mode equations with refinement and norm checks, compare nonzero higher transition terms, integrate switching/radial functions and test zero-momentum/constant-mass controls. Numerical ODE diagnostics supplement the written proof and are not validated integration or formalization. Native, direct science, ordinary and CLI retain unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The free Dirac production report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.165.FREE_FLAT_CLOCK_DIRAC_OUT_ENERGY_AND_COMPLETE_TRANSITION_REMAINDER replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
