"""Read-only replay of the same-source actual smooth weighted extension."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8a_continuation import verify as prior
from p8a_preparation import verify as original

from . import bounds, independent, regularity, weighted
from . import map as equation

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"weighted-extension.json"
PINS = {"A13": "1569abacc14b4c53c951fc10e07eea7e5f17625812705bfb5709772a26261d97",
        "A11": "8dba8215c5a28e2a6379c09e4d452838449065c396e10888aa733b186799456f"}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serialize(value):
    if isinstance(value, dict):
        return {name: serialize(item) for name, item in value.items()}
    if isinstance(value, bool):
        return value
    return str(value)


@cache
def prior_checks():
    if sha(prior.REPORT) != PINS["A13"] or sha(original.REPORT) != PINS["A11"]:
        raise ValueError("a pinned actual-map or causal-inverse input changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return PINS.copy()


def exact_checks():
    groups = {"full_forced_map": equation.identities(), "weighted_prefix_and_gains": weighted.identities(),
              "same_flat_start_regularity": regularity.identities()}
    for name, group in groups.items():
        if any(sp.simplify(value) != 0 for value in group.values()):
            raise ValueError("a weighted-extension exact identity failed: "+name)
    return {name: dict.fromkeys(group, "0") for name, group in groups.items()}


def checked_constants():
    actual = serialize(bounds.calibration())
    separate = independent.replay(json.loads(original.REPORT.read_text()), json.loads(prior.REPORT.read_text()))
    if actual != separate:
        raise ValueError("independent Fraction extension constants disagree")
    return {"complete_weighted_map": actual, "independent_Fraction_replay": separate,
            "smoothness": serialize(regularity.calibration())}


def build_report():
    inputs = prior_checks()
    identities, constants = exact_checks(), checked_constants()
    paths = sorted(ROOT.glob("src/p8a_extension/*.py"))+sorted(ROOT.glob("tests/*.py"))
    paths += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1, "claim": "P8-A.14", "date": "2026-09-06",
        "status": "ACTUAL_SAME_SOURCE_SMOOTH_SEE_EXTENSION_CERTIFIED; FRESH_QSEI_FOCUSING_AND_P8_OPEN",
        "prior_sha256": inputs, "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in paths},
        "formulation": "FORMULATION.md", "full_contraction": "notes/contraction.md",
        "smooth_state_constraint_and_uniqueness": "notes/regularity.md", "source_audit": "notes/sources.md",
        "exact_residuals": identities, "derived_constants": constants,
        "theorem": {
            "field_state_scheme": "same massless minimal scalar, original transported vacuum, fixed radiation, lambda=2sqrt(2) A eta_star^2, gamma=0, physical epsilon=1, delta=10^-14",
            "original_source_window": "L0=10^-10; the identical chosen cutoff chi_L0 and fixed old density defect cbar",
            "new_future_length": "L=10^-6, a factor10000 increase without extending source support",
            "actual_source_free_domain": "x0+L0/2<x<x0+L",
            "same_initial_metric_state_and_constraint": True,
            "positive_Hadamard_state_transported_on_new_metric": True,
            "actual_full_density_and_trace_SEE": True,
            "smooth_on_the_same_longer_slab": True,
            "weighted_ball": "sup exp(-2e6 t)|u'-ubar'|<=2e-7, X(0)=0",
            "whole_history_caps": "duration<=3, |u'|<=10^-5; no state reset at x0",
            "pointwise_geometry": "2<=a<=3,1/3<=h<=1/2, |u|<3e-12",
            "actual_auxiliary_caps": "|a^2 q'|<10^-5, |q|<2e-8, independently established before use",
            "full_pre_inverse_Lipschitz": "less than9e-6, with actual mode/source/anomaly/rolling/local terms",
            "weighted_inverse": "240/769, all A13 causal poles retained",
            "contraction": "less than3e-6",
            "fixed_point_bounds": "weighted<9e-8; pointwise<81e-8 using exp(2)<9",
            "higher_jet_block": "less than10^-4 uniformly in derivative order; no uniform all-jet numerical cap claimed",
            "agrees_with_A11_on_original_slab": True,
            "agreement_reason": "identical original map/data/source and inclusion in its old uniqueness ball",
        },
        "negative_controls": {
            "cutoff_rescaled_to_L": False,
            "weighted_ball_substituted_for_pointwise_barrier": False,
            "full_Einstein_constant_omitted": False,
            "source_pressure_replaced_by_radiation_pressure": False,
            "new_no_runaway_or_future_boundary_condition": False,
            "old_QSEI_or_positive_reference_credit_automatically_transferred": False,
        },
        "verification_boundary": [
            "A13 and A11 hashes checked, A13 and its full A11-to-A1 lineage replayed read-only",
            "Every weighted map constant and pointwise barrier independently reconstructed with Fraction",
            "Independent covariant/source/prefix and false-extrapolation controls included",
            "Infinite-frequency response, Banach, smoothness and actual Hadamard/state/constraint arguments are written proofs, not finite sampling or formalized PDE theorems",
        ],
        "not_established": [
            "a fresh all-sampler QSEI or positive reference credit on the extended domain",
            "macroscopic continuation, runaway suppression, long-time shadowing or optimal domain length",
            "physical validity or small stress fluctuations during the original preparation switch",
            "realistic massive/interacting fields or cosmological incompleteness",
            "completion of P8a or P8",
        ],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("the weighted actual-extension certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8(a) A.14: same-source actual smooth SEE extension replay passed; new QSEI/focusing and P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
