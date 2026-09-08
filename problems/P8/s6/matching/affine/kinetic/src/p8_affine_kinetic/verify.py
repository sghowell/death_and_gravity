"""Read-only specified-affine-kinetic-action certificate replay."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as parent

from . import bridges, curvature, scalar, vector

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"kinetic-promotions.json"
PARENT_SHA = "aeb1f9c7c58fd1dc2a641159a820f3c3fc30265d25eb527530203c2080cec2f4"
sha, serialize = parent.sha, parent.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_affine_kinetic/*.py"))
            +sorted(ROOT.glob("tests/*.py"))+sorted(ROOT.glob("*.md"))
            +sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen exact affine parent report changed")
    expected = json.loads(parent.REPORT.read_text())
    parent.validate_report(expected, parent.build_report())
    return {"S6_37_fully_rebuilt": PARENT_SHA,
            "original_target_and_adopted_contract_rechecked": True,
            "original_physical_metric_and_free_chi_unchanged": True,
            "no_frozen_ancestor_modified": True}


def controls():
    wrong = (0, True, False, 1.0, sp.Float(1), "1", sp.I, sp.oo, -sp.oo, sp.nan)
    calls = [lambda value=value: curvature.sign_control(value) for value in wrong]
    calls += [lambda value=value: vector.require_domain(value)
              for value in (0, -1, True, 1.0, "1/2", sp.I, sp.oo,
                            sp.sqrt(2)/4, sp.Rational(47, 100), sp.Rational(11, 20))]
    calls += [lambda value=value: serialize(value)
              for value in (1.0, sp.Float(1), sp.oo, sp.nan)]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("A forbidden rank/domain/inexact certificate input was accepted")
    if not curvature.sign_control(1)["one_strictly_negative"]:
        raise ValueError("The positive-sign trace-form negative control failed")
    if not curvature.sign_control(-1)["one_strictly_negative"]:
        raise ValueError("The negative-sign trace-form negative control failed")
    return {"rejected_inputs": rejected, "both_trace_form_signs_tested": True,
            "auxiliary_mass_sign_not_a_ghost_test": True,
            "physical_punctured_configuration_test_not_center_inertia_only": True,
            "zero_kinetic_coupling_not_continued_across_rank_change": True,
            "momentum_threshold_not_called_an_EFT_frequency_cutoff": True,
            "three_named_actions_not_a_general_metric_affine_no_go": True}


@cache
def build_report():
    pins = prior_checks()
    groups = {"literal_curvature_screens": curvature.checks(),
              "literal_vector_Schur_and_actual_background": vector.checks(),
              "independent_coupled_scalar_constraints": scalar.checks(),
              "independent_premise_and_unit_interfaces": bridges.checks()}
    residuals = {name: parent.certify_residuals(values) for name, values in groups.items()}
    proofs = {"scalar_domains_and_inertia": scalar.proof_checks(),
              "whole_tube_vector_and_units": bridges.proof_checks()}
    if not all(value is True for group in proofs.values() for value in group.values()):
        raise ValueError("A domain or continuous inequality proof check failed")
    return {
        "schema": 1, "claim": "P8-S6.38.AFFINE_KINETIC", "date": "2026-09-07",
        "status": "THREE_SPECIFIED_KINETIC_PROMOTIONS_SCREENED; GENERAL_AFFINE_UV_AND_ORIGINAL_P8_OPEN",
        "prior_sha256": pins,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": ["notes/curvature.md", "notes/vector.md", "notes/scalar.md", "notes/interfaces.md"],
        "exact_residuals": residuals,
        "named_exact_check_count": sum(map(len, groups.values())),
        "checked_scalar_entries": sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1
                                      for group in groups.values() for value in group.values()),
        "proof_checks": proofs,
        "homothetic_verdict": "Exactly CD/M1 plus a gravitationally coupled Maxwell spectator; the 60-dimensional quotient remains auxiliary",
        "trace_form_verdict": "Opposite signs of unconstrained homogeneous STF-spin3 and axial kinetic terms reject this one literal matrix-trace curvature-square action for either nonzero overall sign",
        "quotient_vector": serialize(vector.calibration()),
        "coupled_scalar": serialize(scalar.calibration()),
        "nonunit_and_parent_interfaces": serialize(bridges.calibration()),
        "controls": controls(),
        "not_established": ["A no-go for arbitrary affine kinetic invariants or their tuned combinations",
                            "An unhealthy low-frequency EFT domain or quantitative ghost-frequency cutoff",
                            "A new vacuum extension, loop or threshold bound, Regge/IR estimate or V/G/B closure",
                            "UV completion, original P8 closure or a new prerequisite for the scoped photon objective"],
        "verification_boundary": "Exact symbolic action/constraint algebra with written whole-domain, representation and physical configuration-velocity proofs; not proof-assistant formalization",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The specified affine kinetic report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.38.AFFINE_KINETIC replay passed; named-action screens only, original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
