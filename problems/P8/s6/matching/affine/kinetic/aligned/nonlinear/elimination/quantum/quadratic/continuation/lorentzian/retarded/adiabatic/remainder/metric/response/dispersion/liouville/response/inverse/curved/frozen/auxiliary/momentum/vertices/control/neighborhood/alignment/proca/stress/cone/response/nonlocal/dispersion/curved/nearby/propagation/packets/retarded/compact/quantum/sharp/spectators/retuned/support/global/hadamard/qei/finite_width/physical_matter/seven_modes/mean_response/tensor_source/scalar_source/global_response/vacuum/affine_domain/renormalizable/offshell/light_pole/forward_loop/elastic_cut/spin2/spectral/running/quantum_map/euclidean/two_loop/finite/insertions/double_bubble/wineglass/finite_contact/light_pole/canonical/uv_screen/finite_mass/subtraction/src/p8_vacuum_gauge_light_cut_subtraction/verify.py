"""Read-only finite-cut subtraction and coefficient enclosure checkpoint."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_finite_mass_gauge_cut import verify as parent

from . import analytic, audit, calibration, domain, logarithm, subtraction

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-first-gauge-light-cut-subtraction.json"
)
PARENT_SHA = "1609098903d12395542aea9b57a2c0ac92e26ddc14746e69aef833d24622b6fe"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_gauge_light_cut_subtraction/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen finite-mass first gauge cut changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_129_fully_rebuilt": PARENT_SHA,
        "same_GY14_unbroken_candidate_and_first_two_gauge_channel": True,
        "old_scalar_action_and_reports_unchanged": True,
        "explicit_low_cut_subtraction_not_an_action_counterterm_retuning": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A finite light-cut subtraction gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.130.FIRST_TWO_GAUGE_LIGHT_CUT_SUBTRACTION",
        "date": "2026-09-10",
        "status": "EXPLICIT_FIRST_TWO_GAUGE_LOW_CUT_REMOVAL_AND_FINITE_MASS_SUBTRACTION_COEFFICIENT_ENCLOSURE; NOT_FULL_MATCHED_AMPLITUDE_ALL_LOOP_ANALYTICITY_POSITIVITY_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/domain.md",
            "notes/subtraction.md",
            "notes/bounds.md",
            "notes/ambiguity.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": len(rows),
        "proof_checks": gates,
        "complex_finite_mass_transition_domain": serialize(payload(domain.data())),
        "finite_cut_and_boundary_germs": serialize(payload(subtraction.data())),
        "holomorphic_coefficient_error_bound": serialize(
            {
                "analytic": payload(analytic.data()),
                "logarithm": payload(logarithm.data()),
            }
        ),
        "same_candidate_subtraction_coefficient": serialize(
            payload(calibration.data())
        ),
        "controls": serialize(audit.controls()),
        "verdict": "For the unchanged prospective GY14-unbroken candidate and its first two-gauge cut at three parent loops, the entire finite-mass transition tail is bounded on a complex s disc. The explicit Cauchy integral over cut invariant zero to six retains both crossed boundary signs and removes exactly this known low-cut discontinuity locally. Its divided-difference representation has controlled analytic boundary germs, so a Cauchy estimate bounds its center second coefficient without differentiating a singular principal-value integral naively. The leading coefficient is K(2 log(2)-21/4), and the complete finite-mass correction is below 10^-190 K at the named mass. Thus the subtraction coefficient lies strictly between -4K and -3K, with its absolute value below 4 times 10^-1620. The regular matched amplitude and high-energy contour remain separate: identical cut data do not fix them. No full amplitude positivity, all-loop mass-gap restoration, finite-gravity or common-parent verdict is claimed.",
        "not_established": [
            "The remaining regular part of the new-model matched amplitude and its b2",
            "All higher-loop light cuts or nonperturbative gauge spectrum",
            "A complete new-model canonical pole or forward-amplitude error budget",
            "A justified high-energy dispersion contour and strict positivity verdict",
            "Finite-gravity IR/Regge Delta and the controlled common bounce parent",
            "Original V G B or P8 closure or exclusion of any full ladder row",
        ],
        "verification_boundary": "Exact inherited normalization and routing constants, explicit Cauchy primitives and crossed germs, removable polynomial quotient diagnostics, rational log enclosures and coefficient-error calibration. Uniform complex-domain, symmetry descent, Ward soft endpoint, Plemelj/Morera gluing and Cauchy arguments are written proofs, not proof-assistant formalized or nonperturbative. Native and ordinary scientific replays use unmodified SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The first gauge light-cut subtraction report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.130.FIRST_TWO_GAUGE_LIGHT_CUT_SUBTRACTION replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
