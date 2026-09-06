"""Read-only replay of the actual short-slab conserved-preparation theorem."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8a_existence import verify as a10
from p8a_remainder import verify as a7
from p8a_residual import verify as a8

from . import actual_map, bounds, independent, regularity

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"smooth-preparation.json"
PINS = {"A10": "c37ee6cb74eb8383abbddf2f240c28676671e57b9f114dd8ec2195267de2d3ce",
        "A8": "dbd88193542224b3eeb76e18923744f3317791ad46961caca3562a1525ba56d1",
        "A7": "cc4ec02bcb23e0ef01d0fb58ad4b95a3bb602a5bbbfc6871d044b2a09a78e963"}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serialize(value):
    if isinstance(value, dict):
        return {name: serialize(item) for name, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(item) for item in value]
    if isinstance(value, bool):
        return value
    return str(value)


@cache
def prior_checks():
    for label, module in (("A10", a10), ("A8", a8), ("A7", a7)):
        if sha(module.REPORT) != PINS[label]:
            raise ValueError("a pinned preparation input changed: "+label)
    a10.validate_report(json.loads(a10.REPORT.read_text()), a10.build_report())
    return PINS.copy()


def exact_checks():
    groups = {"full_actual_trace_constraint_and_conserved_source": actual_map.identities(),
              "complete_shared_history_map": actual_map.difference_identities(),
              "smoothness_majorant_allocation": regularity.identities()}
    for name, group in groups.items():
        if any(sp.simplify(value) != 0 for value in group.values()):
            raise ValueError("an exact preparation identity failed: "+name)
    return {name: dict.fromkeys(group, "0") for name, group in groups.items()}


def checked_constants():
    data = bounds.calibration()
    smooth = regularity.highest_jet_bound(bounds.HISTORY_DERIVATIVE_CAP, bounds.HISTORY,
                                          bounds.LOCAL_CAP, bounds.INVERSE_CAP)
    fractions = independent.replay(json.loads(a7.REPORT.read_text()), json.loads(a8.REPORT.read_text()))
    bridge = {**data["A8_inputs"],
              "center_P": data["center"]["P_difference"],
              "center_q": data["center"]["q_difference"], "center_rhs": data["center"]["rhs"],
              **{key: data["pairs"][key] for key in
                 ("Einstein_P", "curvature_P", "source_P", "auxiliary_Sprime", "local_Wick_terms")},
              "pair_P": data["pairs"]["P"], "pair_q": data["pairs"]["q"],
              "response_quadratic": data["mode_response"]["quadratic"],
              "response_higher": data["mode_response"]["higher"],
              "rhs_lipschitz": data["rhs_lipschitz"], "inverse_bound": data["inverse_bound"],
              "actual_contraction": data["actual_gate"]["contraction"],
              "actual_center_image": data["actual_gate"]["center_image"],
              "actual_self_map": data["actual_gate"]["self_map"],
              "actual_distance": data["actual_gate"]["fixed_point_distance"],
              "auxiliary_P_derived": data["auxiliary_P_derived"],
              "auxiliary_q_derived": data["auxiliary_q_derived"],
              "rounded_self_map": data["rounded"]["self_map"],
              "highest_jet_contraction": smooth["highest_jet_contraction"]}
    bridge.pop("U0_per_delta")
    bridge.pop("U1_per_delta")
    if serialize(bridge) != fractions["constants"]:
        raise ValueError("independent Fraction integration disagrees with the full-map constants")
    if smooth["highest_jet_contraction"] >= bounds.CONTRACTION_CAP:
        raise ValueError("the common all-order smoothness block exceeds its stated cap")
    return {"full_map": serialize(data), "smooth_highest_block": serialize(smooth),
            "independent_Fraction_replay": fractions}


def build_report():
    inputs = prior_checks()
    identities, constants = exact_checks(), checked_constants()
    paths = sorted(ROOT.glob("src/p8a_preparation/*.py"))+sorted(ROOT.glob("tests/*.py"))
    paths += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1, "claim": "P8-A.11", "date": "2026-09-06",
        "status": "SMOOTH_ACTUAL_LOCAL_SEE_WITH_CONSERVED_JOINT_PREPARATION_CERTIFIED; WHOLE_TARGET_AND_P8_OPEN",
        "prior_sha256": inputs,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in paths},
        "formulation": "FORMULATION.md",
        "written_proofs": ["notes/construction.md", "notes/contraction.md", "notes/regularity.md"],
        "source_audit": "notes/sources.md", "remaining_obligations": "notes/remaining.md",
        "exact_residuals": identities, "derived_constants": constants,
        "theorem": {
            "physical_field": "one free massless minimally coupled scalar in the original prepared positive quasifree Hadamard vacuum",
            "fixed_prescription": "lambda=2sqrt(2) A eta_star^2, gamma=0; no new Wick shift or explicit curvature-squared gravity",
            "fixed_classical_source": "ordinary radiation rho_rad=3 A^2/(kappa a_phys^4); zero Lambda",
            "physical_amplitude": "epsilon=1, delta=10^-14, kappa hbar=2880pi^2 delta A^2 eta_star^4",
            "common_history": "unchanged A8 through y=5/2 and a following flat-start neighborhood; no final covariance frozen while changing metric jets",
            "external_source": "prescribed c=chi_L cbar; rho_ext=F c/a^4, p_ext=F[c/(3a^4)-c'/(3h a^4)], F=hbar/(pi^2 A^4 eta_star^8)",
            "cutoff": "smooth 0<=chi_L<=1, equals1 on[0,L/4], equals0 on[L/2,L], L=10^-10",
            "source_globally_compact_in_time": False,
            "complete_trace": "q''+2hq'=u/(60delta)-[u^2/4+h^2(u+h^2)/30]/a^2+8c'/(ha^2)",
            "energy_constraint": "C=3a'^2-3-2880delta(P_rho+c); C'=a^4 h Dtrace; C(0)=0 exactly",
            "ball": "X=u'-ubar', X(0)=0, ||X||<=10^-6; all original nonzero mode history retained",
            "full_map_contraction_upper": "3/25", "full_map_center_upper": "63/1000000000",
            "fixed_point_X_upper": "(63/1000000000)/(1-3/25) <72/1000000000",
            "smoothness": "all Picard iterates have common flat start; every highest-jet block bounded by the same strict contraction on the same slab",
            "actual_state": "positive original vacuum transported on the newly constructed smooth metric; Hadamard propagation then identifies the actual mode stress",
            "unforced_SEE_domain": "x0+L/2 < eta/eta_star-1 < x0+L; smooth full pointwise equation and fixed radiation",
            "uniqueness_scope": "the specified data, cutoff and closed ball; not all possible states or preparations",
        },
        "negative_controls": serialize(actual_map.negative_controls()),
        "verification_boundary": [
            "A10/A8/A7 pins and their lineage replayed; no immutable old file is edited",
            "Fraction-only polynomial integration independently checks all new full-map constants",
            "Physical-clock normalization, nonradiation source integration and constraint propagation have separate algebraic controls",
            "Larger intervals and failed self-map assumptions are explicitly rejected by conditional arithmetic tests",
            "Infinite-frequency convergence, smoothness induction, Banach and Hadamard arguments are written analytic proofs, not finite sampled or Lean proofs",
        ],
        "not_established": [
            "whole-A8-target SEE continuation or long-time shadowing/stability",
            "A9 QSEI transfer to the new exact solution or focusing/incompleteness hypotheses",
            "uniform smallness of all higher curvature jets or preparation pressure derivatives",
            "massive/interacting/realistic-field or semiclassical-fluctuation control",
            "completion of original P8a or P8",
        ],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("the smooth-preparation certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8(a) A.11: smooth local actual SEE and conserved preparation replay passed; whole target and P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
