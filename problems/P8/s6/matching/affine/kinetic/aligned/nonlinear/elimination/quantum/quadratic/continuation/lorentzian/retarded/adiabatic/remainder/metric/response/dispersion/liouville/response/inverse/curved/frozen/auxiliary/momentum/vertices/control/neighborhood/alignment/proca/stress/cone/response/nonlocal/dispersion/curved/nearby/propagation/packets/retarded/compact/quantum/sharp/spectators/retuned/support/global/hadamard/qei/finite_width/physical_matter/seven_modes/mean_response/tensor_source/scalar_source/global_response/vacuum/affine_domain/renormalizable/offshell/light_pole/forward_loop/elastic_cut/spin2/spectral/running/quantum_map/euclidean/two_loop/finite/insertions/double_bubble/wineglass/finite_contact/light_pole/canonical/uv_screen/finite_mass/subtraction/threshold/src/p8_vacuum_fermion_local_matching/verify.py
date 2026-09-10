"""Read-only fermion local-reference and finite-functional health checkpoint."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_gauge_light_cut_subtraction import verify as parent

from . import anchors, audit, calibration, potential, reference, twopoint

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-fermion-local-matching.json"
PARENT_SHA = "0698704470509fa022fa4ccf098a27c20b4724b1235eb7387785fc735ce7bb78"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_fermion_local_matching/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen first light-cut subtraction changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_130_fully_rebuilt": PARENT_SHA,
        "same_GY14_unbroken_prospective_reference_boundary": True,
        "old_scalar_action_and_reports_unchanged": True,
        "fermion_local_reference_ledger_not_full_new_model_matching": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A fermion local-reference gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.131.FERMION_LOCAL_REFERENCE_MATCHING",
        "date": "2026-09-10",
        "status": "COMPLETE_ONE_LOOP_FERMION_INCREMENT_AND_SELECTED_TREE_PLUS_FERMION_FUNCTIONAL_REFERENCE_HEALTH; NOT_FULL_NEW_MODEL_QUANTUM_POTENTIAL_AMPLITUDE_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/potential.md",
            "notes/anchors.md",
            "notes/twopoint.md",
            "notes/reference.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": len(rows),
        "proof_checks": gates,
        "all_flavor_fermion_potential": serialize(payload(potential.data())),
        "two_point_and_cancellation_free_anchors": serialize(
            {"two_point": payload(twopoint.data()), "anchors": payload(anchors.data())}
        ),
        "explicit_local_reference_and_field_map": serialize(payload(reference.data())),
        "same_candidate_selected_functional_bounds": serialize(
            payload(calibration.data())
        ),
        "controls": serialize(audit.controls()),
        "verdict": "For the unchanged prospective GY14-unbroken reference boundary, the complete one-loop fermion determinant retains all fourteen vacuum-energy constants, the large finite mass threshold and the quartic threshold. An independent Dirac trace gives the two-point increment and its entire affine dimensional reference. Cancellation-free exact moment series enclose the on-shell anchors, and a complex-domain bound controls the twice-subtracted kernel. Explicit local mass and vacuum-energy references followed by the stated constant field map give a unit-residue mass-one pole and no additional pole in the unit disc for the selected tree-scalar plus one-loop-fermion functional. Its potential is positive away from the origin on the finite fermion-gap field domain: the mass curvature, heavy-eliminated quartic and every higher even fermion coefficient are positive there. The large negative reference mass is not a physical tachyon or a naturalness exclusion. Canonical couplings and potential derivatives are not the unchanged MSbar ray values. Old scalar loops and later mixed or gauge loops are not included, so no complete new-model quantum-vacuum, forward-amplitude, V/G/B or original P8 closure verdict follows.",
        "not_established": [
            "The complete new-model quantum potential and all-momentum spectrum",
            "The full four-scalar amplitude and the new-model higher-loop reference forest",
            "Bounds for all omitted terms needed for a strict dispersion verdict",
            "High-energy contours and nonperturbative gauge spectrum",
            "Finite-gravity IR/Regge Delta and the common bounce field/state dictionary",
            "Original V G B or P8 closure or exclusion of any full ladder row",
        ],
        "verification_boundary": "Exact determinant derivatives, Dirac trace, moment recurrences, grouped anchor coefficients, on-shell subtraction identities, field maps and rational calibration. Complete series tails, dimensional reference organization, complex holomorphy, pole exclusion and finite-field positivity are written proofs, not proof-assistant formalized, nonperturbative or full new-model claims. Native and ordinary scientific replays use unmodified SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The fermion local-reference report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.131.FERMION_LOCAL_REFERENCE_MATCHING replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
