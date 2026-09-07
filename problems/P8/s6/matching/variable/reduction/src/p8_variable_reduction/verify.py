"""Read-only formal own-f/CD mismatch certificate with recursively rebuilt ancestry."""
import argparse
import hashlib
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

import sympy as sp
from p8_variable_beta import verify as prior

from . import bridges, dictionary, stationary, tensor

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"stationary-cd-matching.json"
P8 = next(path for path in ROOT.parents if path.name == "P8")
PARENT_SHA = "335cd52028baf56b30c75377aa7e6edd7ec86db2799f49b13f467242674fc4d0"
CONTRACT_SHA = "d00df35d5821b05baa61e897725a22ed1ac21d0536b74bdaece3ab0093215901"
FORMULATION_SHA = "73dd18f7a4b56a0205f9fbc6a4b09b213a3c280b3f0274414a255635b04e2b09"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serialize(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, sp.Basic):
        if value.has(sp.Float, sp.oo, sp.zoo, sp.nan, sp.I):
            raise TypeError("Only exact real formal report expressions are admitted")
        return str(value)
    if isinstance(value, dict):
        return {str(key): serialize(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(item) for item in value]
    if value is None or isinstance(value, (str, int, bool)):
        return value
    raise TypeError("Only exact formal data or JSON metadata are admitted")


def source_files():
    return (sorted(ROOT.glob("src/p8_variable_reduction/*.py"))
            +sorted(ROOT.glob("tests/*.py"))+sorted(ROOT.glob("*.md"))
            +sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    pins = {"S6_20_actual_parent_recursively_rebuilt": (prior.REPORT, PARENT_SHA),
            "adopted_V_G_B_matching_contract": (P8/"s6"/"FORMULATION.md", CONTRACT_SHA),
            "original_physical_frame_and_operator_basis": (P8/"FORMULATION.md", FORMULATION_SHA)}
    if any(sha(path) != digest for path, digest in pins.values()):
        raise ValueError("A frozen parent, contract or physical-frame formulation changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    target = bridges.target_input()
    return {**{key: digest for key, (_, digest) in pins.items()},
            "original_CD_covariant_input": target["sha256"],
            "original_CD_source_hashes_checked": target["source_hash_count"]}


def controls():
    calls = [lambda value=value: dictionary.center(x=value) for value in
             (True, False, 1.0, sp.Float("1"), sp.oo, -sp.oo, sp.nan, "1", sp.Symbol("x"))]
    calls += [lambda value=value: dictionary.center(c=value) for value in (1, 2, 5)]
    calls += [lambda value=value: dictionary.center(x=value)
              for value in (sp.Rational(9, 10), sp.Rational(11, 10))]
    calls += [lambda: dictionary.center(parent_planck_squared=0),
              lambda: dictionary.center(time_scale=-1),
              lambda: dictionary.remainder_floor(-1), lambda: dictionary.remainder_floor(1),
              lambda: dictionary.invariant_remainder_floor(kinetic_error=1),
              lambda: dictionary.invariant_remainder_floor(invariant_error=-1),
              lambda: dictionary.remainder_floor(target_planck_squared=0),
              lambda: serialize(0.1), lambda: serialize(sp.Float("0.1")),
              lambda: serialize(sp.oo), lambda: serialize(sp.nan), lambda: serialize(sp.I)]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An invalid numeric-domain or report input was admitted")
    return {"rejected_inputs": rejected,
            "formal_symbolic_expressions_admitted_only_as_report_algebra": True,
            "no_actual_c_equals_2_action": True,
            "inverse_X_target_grading_retained": True,
            "background_value_not_substituted_for_C4_or_C1_control": True,
            "curvature_square_and_moving_coefficients_retained": True,
            "physical_g_and_chi_not_reassigned_by_field_redefinition": True,
            "no_gap_cutoff_or_controlled_EFT_claim": True}


@cache
def build_report():
    pins = prior_checks()
    audits = [ROOT/"tests"/"test_reduction_independent_audit.py",
              ROOT/"tests"/"test_reduction_center_independent_audit.py"]
    if any(not path.is_file() for path in audits):
        raise ValueError("Both separately authored physical-tensor and center-jet audits are required")
    groups = {"literal_full_potential_own_f_and_conformal_action": stationary.checks(),
              "scalar_IBP_clock_and_matter_dictionary": dictionary.checks(),
              "literal_physical_tensor_geometry_and_variation": tensor.checks(),
              "actual_parent_and_continuous_CD_bridges": bridges.checks()}
    if any(value != 0 for group in groups.values() for value in group.values()):
        raise ValueError("A formal reduction or target identity failed")
    margins = {**{"stationary_"+key: value for key, value in stationary.domain_margins().items()},
               **{"bridge_"+key: value for key, value in bridges.domain_margins().items()}}
    if any(value <= 0 for value in margins.values()):
        raise ValueError("A continuous algebraic-domain margin failed")
    coefficients = dictionary.coefficients()
    jets = stationary.center_jets()
    polynomial = jets["remainder_u4_minus_10752_cleared"].all_coeffs()
    if any(value < 0 for value in polynomial) or not any(value > 0 for value in polynomial):
        raise ValueError("The continuous fourth-derivative margin polynomial failed")
    return {
        "schema": 1, "claim": "P8-S6.30.REDUCTION", "date": "2026-09-07",
        "status": "FORMAL_SOURCE_PRESERVING_ACTION_AND_NECESSARY_CD_REMAINDER_DEFECT; ORIGINAL_S6_P8_OPEN",
        "prior_sha256": pins,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": ["notes/stationary.md", "notes/dictionary.md", "notes/tensor-obstruction.md"],
        "source_audit": "notes/sources.md",
        "exact_residuals": {name: dict.fromkeys(group, "0") for name, group in groups.items()},
        "exact_residual_count": sum(map(len, groups.values())),
        "strict_domain_constants_and_margins": serialize(margins),
        "stationary_action_calibration": serialize(stationary.calibration()),
        "full_scalar_tensor_normal_form_with_separate_curvature_square": serialize({
            key: value for key, value in coefficients.items() if key != "symbols"}),
        "dictionary_and_necessary_error_calibration": serialize(dictionary.calibration()),
        "invariant_approximate_remainder_floor": serialize(dictionary.invariant_remainder_floor(
            sp.Rational(1, 10), sp.Rational(1, 10))),
        "literal_physical_tensor_obstruction": serialize(tensor.calibration()),
        "physical_g_source_and_matter_contacts": serialize({
            key: value for key, value in dictionary.contacts().items() if key != "symbols"}),
        "actual_background_center_jets": serialize({
            key: value for key, value in jets.items() if key != "remainder_u4_minus_10752_cleared"}),
        "continuous_fourth_jet_margin_polynomial": {
            "definition": "c^2*(R00_uuuu(0,c)-10752), c=2+delta",
            "coefficients_descending_delta": serialize(polynomial),
            "strict_for_all_delta_positive": True,
            "not_a_sampled_time_or_parameter_bound": True},
        "controls": controls(),
        "not_established": [
            "A selected differential inverse, homogeneous state or controlled derivative expansion",
            "A rolling gap, physical cutoff, loop or all-operator omission bound",
            "An action or solution at the excluded c=2 endpoint",
            "A ghost or UV exclusion from the finite higher-derivative truncation",
            "Absence of every resummed, finite-band or different-parent description",
            "Equality of original CD/matter data after silently changing the physical metric",
            "A V/G pass, whole-row classification with UV, global healthy completion or original P8 closure"],
        "verification_boundary": "Exact literal tensor and potential algebra, actual-input bridges, separately authored scientific audits and written continuous IBP/inequality proof; not proof-assistant formalization",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The stationary-CD matching certificate differs from exact read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.30.REDUCTION: formal own-f/CD mismatch replay passed; original S6/P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
