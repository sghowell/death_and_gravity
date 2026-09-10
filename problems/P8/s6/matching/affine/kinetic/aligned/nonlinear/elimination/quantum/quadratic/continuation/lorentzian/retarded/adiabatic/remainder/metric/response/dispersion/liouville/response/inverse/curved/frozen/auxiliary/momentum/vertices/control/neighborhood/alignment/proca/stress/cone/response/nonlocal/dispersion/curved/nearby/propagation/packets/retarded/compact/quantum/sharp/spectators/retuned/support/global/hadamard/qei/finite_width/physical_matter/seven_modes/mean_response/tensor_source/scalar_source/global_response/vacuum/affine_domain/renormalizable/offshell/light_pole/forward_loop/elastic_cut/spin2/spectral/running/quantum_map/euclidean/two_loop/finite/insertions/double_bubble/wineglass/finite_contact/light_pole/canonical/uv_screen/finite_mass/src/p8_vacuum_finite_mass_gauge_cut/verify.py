"""Read-only finite-mass first two-gauge cut checkpoint."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_gauge_yukawa_screen import verify as parent

from . import audit, calibration, cuts, dirac, kinematics, series

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-finite-mass-first-gauge-cut.json"
PARENT_SHA = "279ec9e7cefbb9b92a3700f92cb3e114ec55c723ee247ecf0b22bbef783027c9"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_finite_mass_gauge_cut/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen leading gauge-Yukawa screen changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_128_fully_rebuilt": PARENT_SHA,
        "same_GY14_unbroken_candidate_and_prospective_boundary": True,
        "old_scalar_action_and_reports_not_modified": True,
        "momentum_remainder_extends_not_relabels_leading_parent_claim": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A finite-mass first gauge-cut gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.129.FINITE_MASS_FIRST_TWO_GAUGE_CUT",
        "date": "2026-09-10",
        "status": "UNIFORM_ALL_DEGREE_ONE_LOOP_FERMION_BOX_MOMENTUM_REMAINDER_AND_NONZERO_COMPLETE_FIRST_TWO_GAUGE_CUT_AT_THREE_PARENT_LOOPS; NOT_ALL_LOOP_AMPLITUDE_CONFINEMENT_UV_NO_GO_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/dirac.md",
            "notes/series.md",
            "notes/kinematics.md",
            "notes/cut.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entries(),
        "proof_checks": gates,
        "Dirac_Ward_and_complex_kinematics": serialize(
            {"Dirac": payload(dirac.data()), "kinematics": payload(kinematics.data())}
        ),
        "all_degree_finite_mass_box_remainder": serialize(payload(series.data())),
        "first_two_gauge_cut_bound": serialize(payload(cuts.data())),
        "same_candidate_rational_calibration": serialize(payload(calibration.data())),
        "controls": serialize(audit.controls()),
        "verdict": "For the same prospective GY14-unbroken candidate, all Taylor degrees at least three of each one-loop fermion transition box are uniformly bounded after the full regulated Ward cancellations. The independent rank-14 Ward system fixes the degree-two tensor and the inherited determinant fixes its coefficient. Sewing the two holomorphically continued transitions gives a complete first two-gauge cut bound, with opposite crossed-channel boundary signs. At s=3/2 and mF=10^200, the relative finite-mass error is strictly below 5 times 10^-194, so the boundary mismatch inside the old massive scalar disc cannot vanish at this perturbative order. The crossing center and inconclusive bounds are not misclassified. Direct inheritance of the old gapped dispersion hypotheses remains unjustified. No all-loop amplitude, confinement spectrum, entire candidate or row exclusion, finite-gravity or common-parent result is asserted.",
        "not_established": [
            "All higher-loop or nonperturbative corrections to the gauge candidate",
            "A confinement mass gap or exact gauge-sector spectrum",
            "A complete new-model three-loop amplitude or b2 error budget",
            "A new applicable infrared-subtracted or gapped dispersion relation",
            "Finite-gravity IR/Regge Delta and controlled common bounce-parent matching",
            "Original V G B or P8 closure or exclusion of any whole ladder row",
        ],
        "verification_boundary": "Exact Clifford and Ward matrices, complete finite-degree tensor-basis linear constraints, six routing inventories, all-degree radial and geometric identities, rational cut normalization and strict error inequalities. Uniform norm estimates, Ward-regulated Taylor subtraction, analytic continuation, channel completeness and cutting arguments are written proofs, not proof-assistant formalized or nonperturbative. Native and ordinary scientific replays use unmodified SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The finite-mass first gauge-cut report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.129.FINITE_MASS_FIRST_TWO_GAUGE_CUT replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
