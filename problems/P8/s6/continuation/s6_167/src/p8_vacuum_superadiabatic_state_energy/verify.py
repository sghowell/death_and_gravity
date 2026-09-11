"""Read-only complete exact-frame transition and local state-energy difference."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_flat_dirac_hadamard import verify as parent

from . import audit, energy, rotation, tube

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-superadiabatic-state-energy.json"
PARENT_SHA = "b7843dda9b7b2a3687b1a5ce8ef2b57eb860f2517ad1a708e3115412a251df82"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_superadiabatic_state_energy/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen Hadamard state parent report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_166_fully_rebuilt": PARENT_SHA,
        "same_free_flat_operator_and_in_out_Hadamard_states": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("An exact-frame state-energy proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.167.EXACT_FOUR_FRAME_DIRAC_TRANSITION_AND_UNIFORM_LOCAL_IN_OUT_ENERGY_DIFFERENCE",
        "date": "2026-09-10",
        "status": "COMPLETE_EXACT_QUADRATIC_TRANSITION_WITH_P_MINUS_FIVE_UV_BOUND_AND_UNIFORM_FREE_FLAT_HADAMARD_IN_OUT_LOCAL_ENERGY_DIFFERENCE; NOT_ABSOLUTE_STRESS_INTERACTING_CURVED_STATE_BACKREACTION_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/complex_disk.md",
            "notes/frames.md",
            "notes/energy.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "analytic_complex_disk_and_exact_frames": serialize(
            {"tube": payload(tube.data()), "rotation": payload(rotation.data())}
        ),
        "complete_transition_and_energy_enclosures": serialize(payload(energy.data())),
        "observable_and_source_context": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "Four exact rotating-frame steps preserve the asymptotic state and retain the final off-diagonal connection. A fixed complex disk, strict frequency gaps and Cauchy estimates bound the exact transition by 2^42 pDelta/(tau^4 E^6), not by a truncated adiabatic evolution. The actual uniform amplitude is below 1e-390 and out-particle energy below 1e20. The same-operator Hadamard in/out local T00 difference, including coherence, is uniformly below 1e411 at every real time. These are quadratic flat-space statements, not absolute renormalized stress, a physical interacting cutoff or bounce backreaction.",
        "not_established": [
            "An absolute local renormalized stress bound or a relative error against bounce density",
            "An interacting gauge state, curved FLRW state, backreaction or common parent background solution",
            "The physical EFT cutoff, omitted quantum loops, quantum target matching or global V/G/B",
            "Original P8 closure or formalization of the written analytic arguments",
        ],
        "verification_boundary": "Exact frame/connection and projector algebra, a written fixed-disk induction and unitarity bound, complete momentum integration, and exact actual-parameter inequalities. Independent Taylor-jet and direct derivative diagnostics test the full frame recurrence at 1200-digit precision for the actual parameters; projector coherence, the exact mode energy-exchange identity and radial quadratures test normalization and scope. Numerical diagnostics are not validated integration. The external references supply context and the inherited Hadamard theorem, not these new constants. Native, direct science, ordinary and CLI use unmodified SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The exact-frame state-energy report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.167.EXACT_FOUR_FRAME_DIRAC_TRANSITION_AND_UNIFORM_LOCAL_IN_OUT_ENERGY_DIFFERENCE replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
