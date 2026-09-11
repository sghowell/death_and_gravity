"""Read-only actual curved free Dirac Hadamard state and state-difference report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_flat_dirac_pressure import verify as parent

from . import audit, energy, geometry, hadamard, transition, tube

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-curved-dirac-state.json"
PARENT_SHA = "97f1821d03f6e2f60367caf25bbc650b7122e31aa07cb77409755c084b89a972"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_curved_dirac_state/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen homogeneous flat pressure report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_169_fully_rebuilt": PARENT_SHA,
        "same_mass_profile_and_original_physical_frame": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("An actual curved Dirac state proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.170.ACTUAL_CD_CURVED_DIRAC_HADAMARD_IN_OUT_STATE_AND_UNIFORM_STATE_STRESS_DIFFERENCE",
        "date": "2026-09-10",
        "status": "CURVED_FREE_HADAMARD_IN_OUT_STATE_AND_UNIFORM_STATE_DIFFERENCE_BOUNDED; NOT_ABSOLUTE_CURVED_INTERACTING_PARENT_STRESS_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/geometry.md",
            "notes/disk.md",
            "notes/frames.md",
            "notes/hadamard.md",
            "notes/energy.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "actual_geometry_exact_frames_and_Hadamard_state": serialize(
            {
                "geometry": payload(geometry.data()),
                "disk": payload(tube.data()),
                "transition": payload(transition.data()),
                "Hadamard": payload(hadamard.data()),
            }
        ),
        "uniform_physical_state_stress_difference": serialize(payload(energy.data())),
        "observable_and_source_context": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The prescribed actual CD geometry has exact free in/out states for all42 fermion copies. Variable disks bound arbitrary finite exact rotations; all-order comparison with the checked local external reference identifies these SAME states as Hadamard. Twenty frames give uniform beta below1e-1890 and local energy-density difference below1e-1090, with pressure difference at most one sixth this conservative bound. These are state differences, not absolute curved or full parent stress, a controlled parent bounce or original P8 closure.",
        "not_established": [
            "Absolute curved stress, physical Newton/curvature references, vacuum polarization or anomaly bounds",
            "Full canonical scalar/heavy/gauge interacting parent state, physical cutoff and loop remainders",
            "Quantum target matching, controlled parent background/response, V/G/B or original P8 closure",
            "Momentum-uniform Moller convergence, global Fock implementability or formalization of external theorems",
        ],
        "verification_boundary": "Written variable-disk, all-order exact-frame, tail and Hadamard-comparison proofs; exact geometry, Pauli/Sylvester and full radial algebra; independent full-function jets, complex disks, high-precision actual-scale refinements and complete numerical quadratures. Finite tests are not validated integration or a formal proof. Native, direct science, ordinary and CLI use original SymPy; only full regression uses its separately audited exact GCD adapter.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The actual curved Dirac state report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.170.ACTUAL_CD_CURVED_DIRAC_HADAMARD_IN_OUT_STATE_AND_UNIFORM_STATE_STRESS_DIFFERENCE replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
