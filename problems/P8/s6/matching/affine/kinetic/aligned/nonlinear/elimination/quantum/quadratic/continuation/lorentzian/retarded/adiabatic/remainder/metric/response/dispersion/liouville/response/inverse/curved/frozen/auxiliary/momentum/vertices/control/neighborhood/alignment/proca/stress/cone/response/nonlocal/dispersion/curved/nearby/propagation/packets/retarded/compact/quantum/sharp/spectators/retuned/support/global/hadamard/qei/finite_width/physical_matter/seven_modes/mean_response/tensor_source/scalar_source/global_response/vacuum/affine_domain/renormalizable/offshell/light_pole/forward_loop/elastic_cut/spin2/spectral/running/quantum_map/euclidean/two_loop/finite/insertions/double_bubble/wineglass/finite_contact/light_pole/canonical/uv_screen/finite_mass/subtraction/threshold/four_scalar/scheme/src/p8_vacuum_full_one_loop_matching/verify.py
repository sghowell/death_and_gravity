"""Read-only full scalar-fermion one-loop canonical reference matching."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_fermion_four_scalar import verify as parent

from . import amplitude, audit, calibration, local, logarithm, scheme

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-full-one-loop-matching.json"
PARENT_SHA = "abbc5dfd0adffad21deb187280c1da862845c838ab023c87c8713166320f935d"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_full_one_loop_matching/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen fermion four-scalar report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_132_fully_rebuilt": PARENT_SHA,
        "same_prospective_GY14_reference_boundary": True,
        "old_scalar_action_and_reports_unchanged": True,
        "complete_one_loop_not_all_loop_error_or_P8_closure": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A complete one-loop matching gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.133.FULL_ONE_LOOP_CANONICAL_MATCHING",
        "date": "2026-09-10",
        "status": "COMPLETE_SCALAR_FERMION_ONE_LOOP_CANONICAL_REFERENCE_MATCHING; NOT_HIGHER_LOOP_ERROR_HIGH_ENERGY_V_FINITE_G_COMMON_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/local.md",
            "notes/scheme.md",
            "notes/amplitude.md",
            "notes/bounds.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": len(rows),
        "proof_checks": gates,
        "full_local_reference_ledger": serialize(payload(local.data())),
        "bare_coupling_and_scheme_conversion": serialize(
            {"scheme": payload(scheme.data()), "logarithm": payload(logarithm.data())}
        ),
        "complete_one_loop_amplitude": serialize(payload(amplitude.data())),
        "actual_pole_and_coefficient_bounds": serialize(payload(calibration.data())),
        "controls": serialize(audit.controls()),
        "verdict": "At the named MS interaction boundary, all one-loop scalar and fermion contributions to the canonical light pole and forward second coefficient are combined in one explicit reference scheme. Independent bare-coupling and total-counterterm identities retain the finite scale logarithm, the previously fixed scalar quartic contact and both sector residue signs. The complete scalar one-loop function is converted from its old canonical subtraction scheme; all first fermion boxes are included without an extra LSZ copy. The large local mass reference, H one-point subtraction and all-species vacuum constant are explicit. Rational logarithm, pole and box enclosures place the complete formal one-loop forward coefficient within a relative 10^-6 error of its positive tree value. The normalized one-loop Phi kernel has a unit-residue mass-one pole, no other pole on the unit disc, and positive local potential curvature. These are full one-loop results, not new-model higher-loop error bounds, global potential stability, exact gauge spectrum, justified V contours, finite-gravity G, common bounce-parent B or original P8 closure.",
        "not_established": [
            "Bounds for the enlarged model's two-loop and all later omitted contributions",
            "Global quantum-potential stability and the complete exact gauge spectrum",
            "High-energy boundedness and a strict full V dispersion verdict",
            "Finite-gravity IR/Regge Delta with justified regulator and contour",
            "A controlled common bounce-parent field/state/cutoff dictionary",
            "Original V G B or P8 closure or exclusion of any full ladder row",
        ],
        "verification_boundary": "Exact bare field/coupling identities, independent local-counterterm conversion, whole one-loop references, complete graph ownership, forward derivatives and rational calibration. Dimensional subtraction, regulator conventions, logarithm tails and the inherited complex pole/amplitude estimates are written proofs, not proof-assistant formalized or all-order QFT claims. Native and ordinary scientific replays use unmodified SymPy with explicit interpreter recursion and trusted-integer formatting allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The full one-loop matching report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.133.FULL_ONE_LOOP_CANONICAL_MATCHING replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
