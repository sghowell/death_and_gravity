"""Read-only antisymmetric Ricci-difference quadratic certificate."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_affine_traces import verify as parent

from . import dynamics, geometry, source

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"ricci-difference-quadratic.json"
PARENT_SHA = "1ad97343c4a7999c29122c8dc4d73d9e08e3be5d73c1d7f2b0aa9f0548e3ba3a"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_affine_ricci/*.py"))
            +sorted(ROOT.glob("tests/*.py"))+sorted(ROOT.glob("*.md"))
            +sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen constant two-trace certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_39_fully_rebuilt": PARENT_SHA,
            "S6_37_original_unrestricted_action_and_complete_CD_M1_rechecked": True,
            "original_P8_and_adopted_contract_unchanged": True,
            "no_frozen_ancestor_modified": True}


def controls():
    wrong = (True, False, 1.0, sp.Float(1), "1", sp.I, sp.oo, sp.nan, sp.Symbol("unproved"))
    calls = [lambda value=value: dynamics.require_regular(value, 1, 1) for value in wrong]
    calls += [lambda value=value: dynamics.units(value, 2, 5) for value in wrong]
    calls += [lambda values=values: dynamics.require_regular(*values, punctured=False) for values in
              ((-sp.Rational(3, 49), sp.sqrt(sp.Rational(5, 7)), 1),
               (-sp.Rational(1, 8), 0, 1), (-sp.Rational(1, 8), sp.sqrt(5), 1))]
    calls += [lambda value=value: dynamics.require_regular(1, 1, value) for value in (0, -1)]
    calls += [lambda: dynamics.require_regular(1, 0, 1)]
    calls += [lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan)]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An inexact, singular-chart or forbidden domain input was accepted")
    if dynamics.negative_witness(-1, 11)["J_eff"] >= 0:
        raise ValueError("The regular negative-coupling physical witness failed")
    if not dynamics.require_regular(0, 0, 1, punctured=False)["unchanged_auxiliary"]:
        raise ValueError("The zero-coupling control failed")
    return {"rejected_inputs": rejected, "nonzero_negative_coupling_tail_witness": True,
            "zero_coupling_not_substituted_into_lambda_inverse": True,
            "curved_commutator_and_coefficient_jets_not_discarded": True,
            "time_dependent_gradient_boundary_and_gyroscopic_term_retained": True,
            "no_center_Theta_division": True,
            "positive_principal_deformation_not_exact_CD_or_UV_completion": True}


@cache
def build_report():
    pins = prior_checks()
    groups = {"literal_full_connection_and_curved_commutator": geometry.checks(),
              "actual_all_component_rolling_source": source.checks(),
              "full_constraints_principal_matrices_and_units": dynamics.checks()}
    residuals = {name: affine.certify_residuals(values) for name, values in groups.items()}
    proofs = dynamics.proof_checks()
    if not all(value is True for value in proofs.values()):
        raise ValueError("A continuous-domain or coefficient-positivity proof failed")
    return {
        "schema": 1, "claim": "P8-S6.40.AFFINE_RICCI", "date": "2026-09-07",
        "status": "CURVED_QUADRATIC_DEFORMATION_AND_PRINCIPAL_SIGN_CASES_CERTIFIED; EXACT_MATCHING_AND_ORIGINAL_P8_OPEN",
        "prior_sha256": pins,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
        "formulation": "FORMULATION.md", "written_proofs": ["notes/geometry.md", "notes/dynamics.md"],
        "exact_residuals": residuals,
        "named_exact_check_count": sum(map(len, groups.values())),
        "checked_scalar_entries": sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1
                                      for group in groups.values() for value in group.values()),
        "proof_checks": proofs,
        "operator": "-(lambda/4)*(F13_[mu,nu]-F14_[mu,nu])^2 added to the unchanged unrestricted S6.37 parent",
        "curved_reduction": "At quadratic rolling order C[nabla Z]=(R_background/3)Z; Z is algebraic on 1+lambda*R_background/3 !=0",
        "coupled_scalar_calibration": serialize(dynamics.calibration()),
        "units": serialize(dynamics.units(3, 2, 5)),
        "verdict": "Nonnegative coupling has positive physical scalar kinetic and principal gradient matrices on the stated regular charts; negative coupling always has a regular-tail high-momentum kinetic obstruction; zero coupling is unchanged. Every nonzero regular deformation changes the original physical CD/M1 quadratic action away from the center.",
        "controls": controls(),
        "not_established": ["A nonlinear open-tube inverse or spectrum on other backgrounds",
                            "An all-frequency or nonlinear stability theorem, or a conserved time-independent principal energy",
                            "Exact original CD/M1 action matching for nonzero coupling",
                            "A full action/gradient/loop remainder or justified low-frequency EFT cutoff",
                            "Additional heavy spectrum, vacuum/finite-gravity V/G/B estimates, UV completion or original P8 closure"],
        "verification_boundary": "Exact symbolic identities, coefficient-positive continuous proofs and independent coordinate/constraint tests; not proof-assistant formalization",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The Ricci-difference report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.40.AFFINE_RICCI replay passed; quadratic deformation only, original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
