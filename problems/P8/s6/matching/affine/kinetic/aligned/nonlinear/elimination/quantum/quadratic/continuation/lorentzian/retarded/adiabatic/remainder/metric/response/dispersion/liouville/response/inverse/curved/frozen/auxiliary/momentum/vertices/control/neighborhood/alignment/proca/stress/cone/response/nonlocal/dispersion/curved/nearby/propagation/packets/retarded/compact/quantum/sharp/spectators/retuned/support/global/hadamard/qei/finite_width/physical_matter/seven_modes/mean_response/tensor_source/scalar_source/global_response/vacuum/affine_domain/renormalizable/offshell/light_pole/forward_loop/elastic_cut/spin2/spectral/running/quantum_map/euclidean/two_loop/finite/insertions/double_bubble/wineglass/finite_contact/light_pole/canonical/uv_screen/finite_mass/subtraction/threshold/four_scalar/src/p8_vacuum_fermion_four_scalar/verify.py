"""Read-only complete one-loop fermion four-scalar increment checkpoint."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_fermion_local_matching import verify as parent

from . import audit, calibration, kinematics, low_degree, reference, series

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-fermion-four-scalar.json"
PARENT_SHA = "2c8ed16ab546b49992e59dc512050e004b28069c3fee9705a831224bfce0d3a7"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_fermion_four_scalar/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen fermion local-reference report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_131_fully_rebuilt": PARENT_SHA,
        "same_prospective_GY14_reference_boundary": True,
        "old_scalar_action_and_reports_unchanged": True,
        "complete_fermion_increment_not_full_new_model_amplitude": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A fermion four-scalar matching gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.132.FERMION_FOUR_SCALAR_INCREMENT",
        "date": "2026-09-10",
        "status": "COMPLETE_ONE_LOOP_FERMION_FOUR_SCALAR_INCREMENT_AND_SELECTED_CANONICAL_COEFFICIENT_BAND; NOT_FULL_NEW_MODEL_SCALAR_SCHEME_MATCHING_HIGHER_LOOP_ERROR_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/low_degree.md",
            "notes/kinematics.md",
            "notes/series.md",
            "notes/reference.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_count(),
        "proof_checks": gates,
        "complete_low_degree_matching": serialize(payload(low_degree.data())),
        "complex_routes_and_all_degree_remainder": serialize(
            {"kinematics": payload(kinematics.data()), "series": payload(series.data())}
        ),
        "explicit_field_and_potential_dictionary": serialize(payload(reference.data())),
        "same_boundary_selected_coefficient_band": serialize(
            payload(calibration.data())
        ),
        "controls": serialize(audit.controls()),
        "verdict": "For the same prospective GY14-unbroken reference boundary, all six cyclic fermion boxes, both orientations and the active color/flavor factor are retained. An independent Dirac trace reproduces the whole regulated zero-momentum quartic and its finite determinant threshold. The complete S4-invariant quadratic momentum term is fixed by the background two-point derivative and is constant on the equal-mass shell. Every higher momentum degree is bounded on a complex forward domain with explicit invariant descent, so Cauchy bounds the box second coefficient. The local potential quartic, residual box and pole field factors are matched without double counting. At the named boundary the box error is less than 10^-604 of the tree coefficient; the formal one-loop fermion increment and the exactly normalized selected coefficient are positive, with the latter above its reference tree value by less than a relative 10^-204. These are selected tree-plus-fermion-functional results, not the full new-model amplitude: old scalar loops, their finite scheme conversion and required later-loop errors remain uncomputed here. No strict V/G/B or original P8 closure follows.",
        "not_established": [
            "The full new-model amplitude with old scalar loops in one matched reference scheme",
            "All later primitive fermion, scalar and gauge loop error bounds",
            "The sign of the fermion box second coefficient alone",
            "High-energy dispersion contours and the exact gauge spectrum",
            "Finite-gravity IR/Regge Delta and the common bounce field/state/cutoff dictionary",
            "Original V G B or P8 closure or exclusion of any full ladder row",
        ],
        "verification_boundary": "Exact independent Dirac trace, complete finite S4 tensor system, background derivatives, labelled routings, all-degree generating function, field-factor identities and rational interval calibration. Regulated subtraction, infinite-series majorization, invariant holomorphy, Cauchy and diagram ownership are written proofs, not proof-assistant formalized or full new-model quantum claims. Native and ordinary scientific replays use unmodified SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The fermion four-scalar report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.132.FERMION_FOUR_SCALAR_INCREMENT replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
