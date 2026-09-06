"""Read-only replay of the new all-Hadamard bound on A14's longer SEE slab."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8a_extension import verify as actual
from p8a_see_qsei import sampling as clock_lemma
from p8a_see_qsei import scattering
from p8a_see_qsei import verify as qsei_prior

from . import focusing, independent, reference, sampling

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"extended-qsei.json"
PINS = {"A14": "93c0ac119240a8c835da52f456d6e7cdc846e872728717de869e8fb823275ec1",
        "A12": "f72698170ee7ca0468b42ceaa098e69e1a26054c6ffdbf6321f1968eaefd6584"}


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
    for name, module in (("A14", actual), ("A12", qsei_prior)):
        if sha(module.REPORT) != PINS[name]:
            raise ValueError("a pinned actual-extension or general QSEI input changed")
        module.validate_report(json.loads(module.REPORT.read_text()), module.build_report())
    return PINS.copy()


def exact_checks():
    groups = {
        "pinned_general_C1_mode_identities": scattering.identities(),
        "pinned_general_two_frequency_moments": clock_lemma.spectral_identities(),
        "pinned_general_actual_proper_clock": clock_lemma.clock_identities(),
        "new_actual_SEE_signed_reference": reference.identities(),
        "new_signed_reference_absorption": sampling.identities(),
        "new_actual_curvature_index": focusing.identities(),
    }
    for name, group in groups.items():
        if any(sp.simplify(value) != 0 for value in group.values()):
            raise ValueError("an extended actual QSEI identity failed: "+name)
    return {name: dict.fromkeys(group, "0") for name, group in groups.items()}


def checked_constants():
    constants = {"fresh_C1_scattering": sampling.scattering_calibration(),
                 "signed_actual_reference": reference.calibration(),
                 "fresh_proper_H2_sampling": sampling.calibration(),
                 "actual_comoving_index": focusing.calibration()}
    constants = serialize(constants)
    separate = independent.replay(json.loads(actual.REPORT.read_text()))
    if constants != separate:
        raise ValueError("independent Fraction replay differs from the new extended QSEI")
    return {**constants, "independent_Fraction_replay": separate}


def build_report():
    inputs = prior_checks()
    identities, constants = exact_checks(), checked_constants()
    paths = sorted(ROOT.glob("src/p8a_extended_qsei/*.py"))+sorted(ROOT.glob("tests/*.py"))
    paths += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1, "claim": "P8-A.15", "date": "2026-09-06",
        "status": "EXTENDED_ACTUAL_SEE_ALL_HADAMARD_H2_QSEI_CERTIFIED; COMOVING_FOCUSING_AND_P8_OPEN",
        "prior_sha256": inputs, "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in paths},
        "formulation": "FORMULATION.md", "written_proof": "notes/proof.md", "sources": "notes/sources.md",
        "exact_residuals": identities, "derived_constants": constants,
        "theorem": {
            "actual_geometry": "the same-source A14 smooth full SEE extension, not a prescribed old background",
            "field_and_reference_state": "one free real massless minimally coupled scalar; original vacuum transported on the actual smooth metric",
            "target_states": "every Hadamard state of the field on that fixed spacetime; no homogeneity, quasifree, zero-mean or target-SEE hypothesis",
            "fixed_prescription": "lambda=2sqrt(2) A eta_star^2, gamma=0, same ordinary radiation and physical epsilon=1, delta=10^-14",
            "physical_scales": "T0=A eta_star^2; kappa hbar=2880pi^2 delta T0^2",
            "domain": "compact support in x0+L0/2<x<x0+L, L0=10^-10, L=10^-6; real proper H2_0 samplers on compact subintervals",
            "original_source_not_rescaled": True,
            "reference_EED_bound": "at least -hbar/(40pi^2 T0^4); no positive reference credit asserted on the extended domain",
            "absolute_QSEI": "integral f^2 E_target d_tau >= -2hbar/(16pi^2) integral |f_proper_second|^2 d_tau",
            "new_reference_penalty_in_derivative_coefficient": "(2/5)L^4; derived afresh from the actual weighted distance",
            "higher_potential_derivative_caps_needed": False,
            "old_short_interval_coefficient_transferred_without_recalibration": False,
        },
        "focusing_boundary": {
            "congruence": "comoving normals to constant cosmic-time surfaces, either orientation, segments within the certified free domain",
            "proper_span_upper": "3e-6 T0",
            "actual_index_lower": "J[g] T0 >=3/tau_hat-tau_hat/8 >3/4 >=-K T0",
            "sufficient_index_trigger_available": False,
            "normal_jacobian_lower": "8/27",
            "sufficient_Q2_over_available_duration_squared_lower": "2/5",
            "other_hypersurface_or_global_no_go": False,
        },
        "verification_boundary": [
            "Pinned A14 and A12 reports and their complete source-hashed ancestor lineages replayed read-only",
            "New reference, scattering, sampler and focusing constants independently reconstructed with Fraction from A14 inputs",
            "Physical Einstein sign, weighted-to-pointwise loss, support, clock, spectral and signed-credit controls independently reviewed",
            "Positive-type, infinite-frequency estimates, density and index implications are written analytic proofs, not numerical frequency sampling or formalized PDE theorems",
        ],
        "not_established": [
            "positive actual reference EED on the extended domain or an optimal QSEI coefficient",
            "long-time or macroscopic SEE continuation, stability, stress fluctuations or preparation EFT validity",
            "cosmological incompleteness, maximal extension, other hypersurfaces or null/boosted geodesics",
            "massive, nonminimal or interacting realistic fields, completion of original P8a or P8",
        ],
    }


def validate_report(expected, current):
    if expected != current:
        raise ValueError("the extended actual-SEE QSEI certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    report = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), report)
        print("P8(a) A.15: extended actual-SEE all-Hadamard H2 QSEI replay passed; focusing and P8 OPEN")
    else:
        print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
