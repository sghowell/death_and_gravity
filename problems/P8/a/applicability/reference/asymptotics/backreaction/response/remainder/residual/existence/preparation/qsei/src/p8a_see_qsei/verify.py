"""Read-only replay of the actual-SEE all-sampler QSEI and focusing boundary."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8a_preparation import verify as prior
from p8a_qsei import verify as positive_type_prior

from . import focusing, independent, reference, sampling, scattering

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"see-qsei.json"
PINS = {"A11": "8dba8215c5a28e2a6379c09e4d452838449065c396e10888aa733b186799456f",
        "A9": "514b4ee50201a13c45dd02ca7bd4879dd58c95453740a0466c5c9c215b04cf20"}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serialize(value):
    if isinstance(value, dict):
        return {key: serialize(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serialize(item) for item in value]
    if isinstance(value, bool):
        return value
    return str(value)


@cache
def prior_checks():
    if sha(prior.REPORT) != PINS["A11"] or sha(positive_type_prior.REPORT) != PINS["A9"]:
        raise ValueError("a pinned actual-solution or positive-type input changed")
    # A.11 replays A.10 -> A.9 -> A.8 -> all earlier A certificates.
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return PINS.copy()


def exact_checks():
    groups = {"exact_scattering_and_endpoints": scattering.identities(),
              "two_frequency_moments": sampling.spectral_identities(),
              "actual_proper_sampler_clock": sampling.clock_identities(),
              "new_SEE_reference_credit": reference.identities(),
              "short_window_index_form": focusing.identities()}
    for name, group in groups.items():
        if any(sp.simplify(value) != 0 for value in group.values()):
            raise ValueError("an actual-SEE QSEI identity failed: "+name)
    return {name: dict.fromkeys(group, "0") for name, group in groups.items()}


def checked_constants():
    modes, samplers = scattering.calibration(), sampling.calibration()
    credit, focus = reference.calibration(), focusing.calibration()
    fractions = independent.replay(json.loads(prior.REPORT.read_text()))
    old = prior.bounds.calibration()
    c0, c1, c2 = modes["forward_two_derivatives"]
    p0, p1 = modes["local_one_derivative"]
    f0, f1, f2 = samplers["sampler_norms"]
    roots = samplers["spectral_error_roots"]
    bridge = {
        "history": modes["history"], "derivative_cap": modes["derivative_cap"],
        "potential_cap": modes["potential_cap"], "length": old["length"], "delta": old["delta"],
        "conformal_span": samplers["conformal_span"], "mode_exponent": modes["modulus_exponent"],
        "q0": modes["q_and_qprime"][0], "q1": modes["q_and_qprime"][1],
        "A0": modes["A_and_two_derivatives"][0], "A1": modes["A_and_two_derivatives"][1],
        "A2": modes["A_and_two_derivatives"][2],
        "forward0": c0, "forward1": c1, "forward2": c2,
        "local0": p0, "local1": p1, "backward_coefficient": modes["backward_coefficient"],
        "bminus1_inverse_k": modes["bminus1_inverse_k"], "bprime_inverse_k": modes["bprime_inverse_k"],
        "qprime_remainder": modes["qprime_minus_uprime_inverse_k"],
        "infrared_error": modes["infrared_error"],
        "sampler0": f0, "sampler1": f1, "sampler2": f2,
        "weighted_sampler": samplers["weighted_first_norm"], "auxiliary_root": samplers["auxiliary_root"],
        "forward_root": roots["forward"], "local_root": roots["local"],
        "history_root": roots["real_history_Parseval"], "remainder_root": roots["history_remainder"],
        "infrared_root": roots["infrared"], "total_root": samplers["total_root"],
        "derived_coefficient": samplers["derived_coefficient"],
        "actual_distance": credit["actual_fixed_point_distance"], "distance_cap": credit["distance_cap"],
        "numerator_loss": credit["numerator_loss_upper"], "positive_EED": credit["positive_EED_dimensionless"],
        "proper_span": focus["proper_span_upper"], "index_lower": focus["index_times_T0_lower"],
        "normal_jacobian": focus["normal_jacobian_lower"],
        "Q2_over_T0_squared": focus["Q2_over_T0_squared"],
        "Q2_over_duration_squared": focus["Q2_over_available_duration_squared_lower"],
    }
    if serialize(bridge) != fractions["constants"]:
        raise ValueError("independent Fraction replay disagrees with the actual-SEE QSEI")
    return {"C1_scattering": serialize(modes), "proper_H2_sampling": serialize(samplers),
            "actual_reference_credit": serialize(credit), "focusing_boundary": serialize(focus),
            "independent_Fraction_replay": fractions}


def build_report():
    inputs = prior_checks()
    identities, constants = exact_checks(), checked_constants()
    paths = sorted(ROOT.glob("src/p8a_see_qsei/*.py"))+sorted(ROOT.glob("tests/*.py"))
    paths += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1, "claim": "P8-A.12", "date": "2026-09-06",
        "status": "ACTUAL_LOCAL_SEE_ALL_HADAMARD_H2_QSEI_CERTIFIED; AVAILABLE_COMOVING_FOCUSING_GATE_FAILS_AND_P8_OPEN",
        "prior_sha256": inputs,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in paths},
        "formulation": "FORMULATION.md", "written_proof": "notes/proof.md",
        "source_audit": "notes/sources.md", "focusing_boundary": "notes/focusing.md",
        "exact_residuals": identities, "derived_constants": constants,
        "theorem": {
            "geometry": "the actual smooth A11 SEE solution, not the old A8 metric",
            "field": "one free real massless minimally coupled scalar",
            "target_states": "all Hadamard states of this field on the fixed new smooth spacetime; no homogeneity, quasifree, zero-mean or target SEE assumption",
            "reference_state": "original positive vacuum transported from the unchanged free past; the actual A11 SEE source",
            "scheme": "lambda=2sqrt(2) A eta_star^2, gamma=0, fixed ordinary radiation; no new finite counterterm",
            "clock_and_amplitude": "x=eta/eta_star-1, T0=A eta_star^2, epsilon=1, delta=10^-14, kappa hbar=2880pi^2 delta T0^2",
            "domain": "sampler support compact in x0+L/2<x<x0+L, L=10^-10; real H2_0 extension on compact subintervals",
            "coefficient": "2 hbar/(16pi^2) multiplying integral |f_proper_second|^2",
            "coefficient_is_globally_optimal": False,
            "positive_reference_EED": "at least hbar/(4199040pi^2 A^4 eta_star^8) on the free half",
            "history_caps": "T<=3, |u'|<=10^-5, |u|<=3*10^-5, actual free-past mode data",
            "higher_potential_jet_norms_assumed": False,
            "new_frequency_argument": "exact forward/backward split, full sampling Parseval, momentum Parseval only for the real u' history, direct infrared Volterra bounds",
            "prior_A9_numerical_bound_transferred": False,
        },
        "focusing_result": {
            "congruence": "comoving normals to constant cosmic-time surfaces, either time orientation, segments inside the certified free half",
            "index_lower": "J[g] T0 >=3/tau_hat-tau_hat/8 >3/4 >=-K T0 for g(0)=1,g(tau)=0",
            "sufficient_index_trigger_available": False,
            "normal_jacobian_lower": "8/27",
            "Q2_over_available_duration_squared": "160000000 for the sufficient coefficient2",
            "global_or_other_hypersurface_exclusion": False,
        },
        "negative_controls": serialize(scattering.controls()),
        "verification_boundary": [
            "A11 and A9 report hashes checked; A11 and its complete lineage replayed read-only",
            "New rational coefficients independently reconstructed from the pinned A11 data with Fraction",
            "Real versus complex Parseval, Born normalization, proper Einstein EED and index-form controls independently audited",
            "Functional analysis, Hadamard positivity, density and index-form implications are written proofs, not finite-frequency sampling or formalized global PDE theorems",
        ],
        "not_established": [
            "QSEI coefficient2 on a longer exact-solution interval or the entire old target",
            "long-time SEE continuation, stability, maximal extension or cosmological incompleteness",
            "other hypersurfaces, boosted/null geodesics, massive/nonminimal/interacting realistic fields",
            "uniform bounds on arbitrary higher curvature jets, fluctuation control or physical EFT validity beyond the model",
            "completion of original P8a or P8",
        ],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("the actual-SEE QSEI certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8(a) A.12: actual-SEE all-Hadamard H2 QSEI replay passed; cosmological focusing and P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
