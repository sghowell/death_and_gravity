"""Read-only regular-root universality and next-order freedom certificate."""
import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_variable_reduction import verify as prior

from . import basis, bridges, general

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"regular-root-freedom.json"
PARENT_SHA = "1095385a6c9da4dded4e6a5234bae202457f9198fd1cd52967c91801906d5977"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_files():
    return (sorted(ROOT.glob("src/p8_general_reduction/*.py"))
            + sorted(ROOT.glob("tests/*.py"))+sorted(ROOT.glob("*.md"))
            + sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(prior.REPORT) != PARENT_SHA:
        raise ValueError("The frozen source-preserving reduction report changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {"S6_30_recursively_rebuilt": PARENT_SHA,
            "original_CD_and_physical_frame": "recursively checked by S6.30",
            "general_interactions_do_not_replace_S6_20_action": True}


def controls():
    beta = (7, 3, sp.Rational(1, 2), sp.Rational(1, 4), -sp.Rational(9, 8))
    bad = (True, False, 1.0, sp.Float("1"), sp.oo, -sp.oo, sp.nan, "2", sp.Symbol("r"))
    calls = [lambda value=value: general.inverse_at(beta, value) for value in bad]
    calls += [lambda: general.inverse_at(beta, 0), lambda: general.inverse_at(beta, -1),
              lambda: general.inverse_at(beta, 2, mg2=0),
              lambda: general.inverse_at(beta, 2, mf2=0),
              lambda: general.inverse_at((0, 1, -sp.Rational(2, 3), sp.Rational(1, 3), 0), 1),
              lambda: general.inverse_at(beta, 1),
              lambda: general.polynomial_at(beta[:4], 2),
              lambda: general.polynomial_at(("7", *beta[1:]), 2)]
    calls += [lambda value=value: prior.serialize(value)
              for value in (0.1, sp.Float("0.1"), sp.oo, sp.nan, sp.I)]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An invalid numeric-domain or report input was admitted")
    return {"rejected_inputs": rejected,
            "symbolic_theorem_not_restricted_to_rational_calibrations": True,
            "P_zero_not_inverted": True, "P_sign_not_a_health_verdict": True,
            "potential_projection_not_full_S6_invariant": True,
            "no_actual_changed_parent_background_or_controlled_EFT_claim": True}


@cache
def build_report():
    pins = prior_checks()
    if not (ROOT/"tests"/"test_general_independent_audit.py").is_file():
        raise ValueError("The independent literal matrix and compact-variation audit is required")
    groups = {"general_potential_stationary_dictionary_and_freedom": general.checks(),
              "explicit_boundary_and_full_Euler_control": basis.checks(),
              "continuous_frozen_parent_specialization": bridges.checks()}
    if any(value != 0 for group in groups.values() for value in group.values()):
        raise ValueError("A general-reduction identity failed")
    table, freedom, bulk = general.curvature_linear_cubic(), general.cubic_freedom(), general.nonboundary_control()
    return {
        "schema": 1, "claim": "P8-S6.31.GENERAL", "date": "2026-09-07",
        "status": "REGULAR_ROOT_RETAINED_ORDER_UNIVERSALITY_AND_GENUINE_NEXT_ORDER_FREEDOM; ORIGINAL_P8_OPEN",
        "prior_sha256": pins,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": ["notes/general-stationary.md", "notes/basis.md"],
        "source_audit": "notes/sources.md",
        "exact_residuals": {name: dict.fromkeys(group, "0") for name, group in groups.items()},
        "exact_residual_count": sum(map(len, groups.values())),
        "calibration": prior.serialize(general.calibration()),
        "full_retained_scalar_normal_form": prior.serialize(general.leading_action()),
        "genuine_formal_S6_deformation": prior.serialize({
            "increments_beta0_to_beta4": freedom["increments"],
            "mixed_h2_cubic_difference": freedom["cubic_shift"],
            "Einstein_second_variation_unchanged": True,
            "full_degree_six_action_individually_computed": False,
            "off_shell_Einstein_delta_L6": bulk["delta_L6"],
            "covariant_metric_Euler_trace": bulk["covariant_metric_Euler_trace"],
            "nonboundary_condition": "eta*Lambda!=0 on the constant-r,P fixture"}),
        "potential_only_curvature_linear_table": prior.serialize({
            "N": table["N"], "coefficients": table["coefficients"],
            "retained_complement": table["retained_complement"],
            "specified_Hessian_free_projection_Xi": table["H_free_projection_Xi"],
            "scope": table["scope"]}),
        "boundary_representative_control": prior.serialize(basis.calibration()),
        "controls": controls(),
        "not_established": [
            "An all-order no-CD theorem or a controlled match to the original CD/M1 action",
            "A selected inverse, homogeneous state, rolling gap, cutoff or all-operator/loop error",
            "Preservation of the full S6.20 rolling solution under the general beta deformation",
            "A global healthy background, finite-gravity positivity/V/G pass or original P8 completion"],
        "verification_boundary": "Written continuous formal-action proof, full symmetric-metric algebra, literal matrix and compact-variation audits; not proof-assistant formalization or a remainder theorem",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The regular-root freedom report differs from exact read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.31.GENERAL replay passed; original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
