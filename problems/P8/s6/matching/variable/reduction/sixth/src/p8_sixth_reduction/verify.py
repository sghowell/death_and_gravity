"""Read-only first-omitted-action tensor and mandatory cancellation replay."""
import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_variable_reduction import dictionary as prior_dictionary
from p8_variable_reduction import verify as prior

from . import bridges, literal, tensor

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"sixth-tensor-cancellation.json"
PARENT_SHA = "1095385a6c9da4dded4e6a5234bae202457f9198fd1cd52967c91801906d5977"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_files():
    return (sorted(ROOT.glob("src/p8_sixth_reduction/*.py"))
            + sorted(ROOT.glob("tests/*.py"))+sorted(ROOT.glob("*.md"))
            + sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(prior.REPORT) != PARENT_SHA:
        raise ValueError("The frozen source-preserving action report changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {"S6_30_recursively_rebuilt": PARENT_SHA,
            "same_S6_20_action_and_original_physical_frame": True,
            "does_not_inherit_a_changed_general_beta_parent": True}


def numeric_center(c=3, parent_planck_squared=1, time_scale=1):
    """Strict exact calibration on the parent's admitted c/M/tau domain."""
    checked = prior_dictionary.center(c=c, parent_planck_squared=parent_planck_squared,
                                      time_scale=time_scale)
    cc, m2, tau = checked["c"], sp.Rational(parent_planck_squared), sp.Rational(time_scale)
    values = tensor.center()
    return {"c": cc, "M_squared": m2, "tau": tau,
            "sixth_E": {n: sp.factor(m2*tau**(n-2)*value.subs(tensor.C, cc))
                        for n, value in values["sixth_E_dimensionless"].items()},
            "fourth_E2": sp.factor(m2*values["fourth_E_dimensionless"][2].subs(tensor.C, cc)),
            "combined_E2": sp.factor(m2*values["combined_E_dimensionless"][2].subs(tensor.C, cc)),
            "probe_is_actual_FLRW_background": False,
            "covariant_Xi_extracted": False}


def controls():
    bad = (True, False, 1.0, sp.Float("1"), sp.oo, -sp.oo, sp.nan, "3", sp.Symbol("c"))
    calls = [lambda value=value: numeric_center(c=value) for value in bad]
    calls += [lambda value=value: numeric_center(c=value) for value in (1, 2, 5)]
    calls += [lambda: numeric_center(parent_planck_squared=0), lambda: numeric_center(time_scale=-1)]
    calls += [lambda value=value: tensor.total_derivative(value)
              for value in (tensor.INV[-1], tensor.LOG[-1], tensor.Q[-1])]
    calls += [lambda value=value: literal._time_derivative(value)
              for value in (literal.A4, literal.L3, literal.Q3, literal.QJ[-1])]
    calls += [lambda value=value: prior.serialize(value)
              for value in (0.1, sp.Float("0.1"), sp.oo, sp.nan, sp.I)]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An invalid numeric/report input or truncated derivative was admitted")
    return {"rejected_inputs": rejected,
            "differentiate_before_center_evaluation": True,
            "fourth_order_cancellation_retained": True,
            "physical_time_not_replaced_by_canonical_scalar": True,
            "unit_volume_probe_not_actual_FLRW": True,
            "individual_floor_not_a_combined_or_all_order_error_bound": True}


@cache
def build_report():
    pins = prior_checks()
    audits = ["test_sixth_independent_audit.py", "test_sixth_center_independent_audit.py"]
    if any(not (ROOT/"tests"/name).is_file() for name in audits):
        raise ValueError("Both separately authored independent scientific audits are required")
    groups = {"full_restricted_action_boundary_Euler_and_cancellation": tensor.checks(),
              "independent_coordinate_ADM_root_series_and_parent_bridges": bridges.checks()}
    if any(value != 0 for group in groups.values() for value in group.values()):
        raise ValueError("An action or actual-profile identity failed")
    data, center = tensor.derive(), tensor.center()
    individual = center["sixth_low_delta_polynomial"]
    combined = center["combined_low_delta_polynomial"]
    if individual.nth(0) != -2 or individual.nth(1) != sp.Rational(7, 4):
        raise ValueError("The continuous individual floor polynomial changed")
    if any(individual.nth(n) >= 0 for n in (2, 3, 4)):
        raise ValueError("The remaining individual polynomial coefficients must be strictly negative")
    if combined.nth(0) != 0 or combined.nth(1) != sp.Rational(7, 2):
        raise ValueError("The required combined center cancellation failed")
    return {
        "schema": 1, "claim": "P8-S6.32.SIXTH", "date": "2026-09-07",
        "status": "FORMAL_FIRST_OMITTED_TENSOR_TERM_AND_MANDATORY_FOURTH_SIXTH_CANCELLATION; ORIGINAL_P8_OPEN",
        "prior_sha256": pins,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
        "formulation": "FORMULATION.md", "written_proofs": ["notes/tensor.md"],
        "source_audit": "notes/sources.md",
        "exact_residuals": {name: dict.fromkeys(group, "0") for name, group in groups.items()},
        "exact_residual_count": sum(map(len, groups.values())),
        "normal_coefficients_over_M_squared": prior.serialize({f"B{n}": data[f"B{n}"] for n in (1, 2, 3)}),
        "full_Euler_coefficients_over_M_squared": prior.serialize(data["euler_coefficients"]),
        "exposed_boundary_APIs": ["derive.covariant_EH_boundary_second", "derive.time_IBP_boundary",
                                  "derive.covariant_to_normal_boundary"],
        "actual_center_B_physical": prior.serialize(center["B_physical_including_M2"]),
        "actual_center_sixth_E_physical": prior.serialize(center["sixth_E_physical"]),
        "actual_center_fourth_E_physical": prior.serialize(center["fourth_E_physical"]),
        "actual_center_combined_E_dimensionless": prior.serialize(center["combined_E_dimensionless"]),
        "continuous_individual_low_symbol_polynomial": prior.serialize(individual.all_coeffs()),
        "mandatory_combined_low_symbol_polynomial": prior.serialize(combined.all_coeffs()),
        "individual_strict_floor": {"domain": "0<delta<=1/100; delta=c-2", "absolute_e2_greater_than": "793/400",
                                     "compact_probe": "q=u^2/2 near zero, smooth compact support", "physical_scale": "M^2/tau^2",
                                     "is_a_combined_or_full_remainder_bound": False},
        "calibration": prior.serialize(tensor.calibration()),
        "exact_numeric_calibration": prior.serialize(numeric_center(4, 3, 2)),
        "controls": controls(),
        "not_established": [
            "A full covariant degree-six Xi or original CD coefficient match",
            "The actual FLRW tensor response, scalar health or a physical cutoff",
            "An order-eight or complete omitted-action error, all-order convergence or exclusion",
            "A specified differential inverse/state/source class, positivity/V/G pass or original P8 closure"],
        "verification_boundary": "Written continuous formal-action and polynomial proof, exact full restricted Euler operator, two separately authored scientific audits; not proof-assistant formalization or a controlled EFT theorem",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The sixth-order tensor cancellation report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.32.SIXTH replay passed; original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
