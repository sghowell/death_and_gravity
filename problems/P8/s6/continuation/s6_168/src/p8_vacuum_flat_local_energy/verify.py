"""Read-only absolute free-flat one-loop local energy with its fixed reference."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_superadiabatic_state_energy import verify as parent

from . import audit, energy, extension, frames, reference, subtraction

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-flat-local-energy.json"
PARENT_SHA = "01ae1949bc7d1b8eb01002e6388c2734d1bfcfb3445cb6bfcb7d7bb96512ddce"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_flat_local_energy/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen exact-frame state-energy report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_167_fully_rebuilt": PARENT_SHA,
        "same_free_flat_operator_exact_frames_and_Hadamard_states": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A fixed-reference absolute local-energy proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.168.FREE_FLAT_DIRAC_ABSOLUTE_ONE_LOOP_LOCAL_ENERGY_WITH_FIXED_SATURATED_MASS_REFERENCE",
        "date": "2026-09-10",
        "status": "BOUNDED_ABSOLUTE_FREE_FLAT_ONE_LOOP_LOCAL_T00_IN_EXPLICIT_MS_VACUUM_AND_SATURATED_MASS_REFERENCE_PRESCRIPTION; NOT_FULL_PHYSICAL_STRESS_CURVED_INTERACTING_STATE_HIGHER_LOOP_REMAINDER_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/frames.md",
            "notes/dimensional.md",
            "notes/reference.md",
            "notes/energy.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "exact_projector_and_dimensional_subtraction": serialize(
            {
                "frames": payload(frames.data()),
                "subtraction": payload(subtraction.data()),
            }
        ),
        "fixed_reference_and_absolute_energy": serialize(
            {
                "reference": payload(reference.data()),
                "extension": payload(extension.data()),
                "energy": payload(energy.data()),
            }
        ),
        "observable_and_source_context": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The full exact quadratic state is compared with four instantaneous rotating frames, not replaced by an adiabatic truncation. Complete momentum estimates bound its subtracted remainder below1e411. Full dimensional MS finite terms and the fixed vacuum reference are restored. The separately named GY14-SAT8-MR finite mass-reference extension preserves all named vacuum coefficients through loop two while bounding the entire paired potential. All local terms and exact state remainders sum to less than1e789, or1e-11 of the named kappa reference, for every real time. This specified free one-loop T00 is not the full physical curved/interacting parent stress or a bounce-density relative error.",
        "not_established": [
            "Full pressure, scalar source, curved Hadamard state or interacting gauge/fermion stress",
            "Complete canonical-field and heavy-background energy or higher-loop and physical-cutoff remainders",
            "Quantum target matching, a common-parent clock or bounce solution, global V contours/cuts or finite-gravity G",
            "Original P8 closure, an identity of full higher-field functionals, or formalization of the analytic proofs",
        ],
        "verification_boundary": "Exact alternating-frame/projector and dimensional gamma finite-part algebra, full radial integrals and a complete convergent paired potential tail. Independent Pauli products, nonzero-regulator radial quadratures, 600-digit actual mass-anchor cancellation, full potential/anchor integrals, curvature/source and unsaturated-counterterm controls supplement the written analytic arguments. Numerical diagnostics are not validated integration. All own exact payloads are serialized before freeze. Native, direct science, ordinary and CLI use unmodified SymPy; only full regression uses the separately audited exact GCD adapter.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The fixed-reference free local-energy report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.168.FREE_FLAT_DIRAC_ABSOLUTE_ONE_LOOP_LOCAL_ENERGY_WITH_FIXED_SATURATED_MASS_REFERENCE replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
