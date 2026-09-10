"""Read-only finite MS quadratic primitive slope certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_fermion_quadratic_forests import verify as parent

from . import audit, calibration, reference, remainder

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-fermion-ms-slopes.json"
PARENT_SHA = "508c59cf4c3a58222bd3afd5af28b242643b6a5fe004e869af335b59e8a4d20f"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_fermion_ms_slopes/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen MS-slope parent report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_145_fully_rebuilt": PARENT_SHA,
        "same_GY14_nongravitational_reference_model": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
        "only_two_primitive_MS_slopes_advance": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("An MS-slope proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.146.FINITE_MS_QUADRATIC_PRIMITIVE_SLOPES",
        "date": "2026-09-10",
        "status": "BOTH_QUADRATIC_PRIMITIVE_NONLOCAL_AND_MS_SLOPE_BOUNDS; NOT_FINITE_MASS_REFERENCE_FULL_MATCHING_CANONICAL_POLE_HIGHER_LOOPS_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/tensors.md",
            "notes/masters.md",
            "notes/forests.md",
            "notes/mass_difference.md",
            "notes/pole.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "exact_massless_MS_reference": serialize(payload(reference.data())),
        "full_scalar_mass_difference": serialize(payload(remainder.data())),
        "actual_primitive_slope_enclosure": serialize(payload(calibration.data())),
        "partial_primitive_frontier": serialize(audit.frontier()),
        "controls": serialize(audit.controls()),
        "verdict": "The complete three-word quadratic primitive tensors are reduced to dimensionally regulated Gaussian masters, with proper fermion mass/kinetic and Yukawa counterterms retained. The massless-boson finite Euclidean slopes are 49/12 for scalar NY^2/Q^2 and -37/6 for gauge NYaC_F/Q^2. Finite epsilon times pole products cancel the pi^2 constants exactly. The actual scalar boson mass is restored by a fully convergent joint integral bound. The prior soft-tail disc bounds the shift from zero momentum to the on-shell slope. Both primitive MS slopes are below 1e-410 in absolute value together. Finite mass references, other matching and canonical contributions, the complete two-loop pole/error and original P8 remain open.",
        "not_established": [
            "Finite MS mass references in the two quadratic primitive rows",
            "Two vacuum primitive rows and other counterterm/matching insertions",
            "Complete matched/canonical two-loop amplitude, pole or error",
            "Higher-loop truncation control, V contours, finite-gravity G, common-parent B or original P8 closure",
        ],
        "verification_boundary": "Exact dimension-symbolic traces, Gaussian-master reduction, independently normalized vacuum mass derivatives, proper-counterterm Laurent expansion, joint finite mass-difference bounds and mutation-checked scope. Numerical checks include explicit matrix derivatives in dimensions three, four and five and nonzero-regulator integrals. Analytic arguments are written proofs, not formalization or independent peer review. Native, ordinary, CLI and direct science retain unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The MS-slope report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.146.FINITE_MS_QUADRATIC_PRIMITIVE_SLOPES replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
