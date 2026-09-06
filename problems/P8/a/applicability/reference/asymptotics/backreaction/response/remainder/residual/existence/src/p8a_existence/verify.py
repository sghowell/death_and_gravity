"""Read-only actual-response continuity and conditional causal-inverse replay."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8a_qsei import verify as prior
from p8a_remainder import verify as remainder

from . import independent, log_inverse, mode_lipschitz

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"response-inverse.json"
PRIOR_SHA = "514b4ee50201a13c45dd02ca7bd4879dd58c95453740a0466c5c9c215b04cf20"
A7_SHA = "cc4ec02bcb23e0ef01d0fb58ad4b95a3bb602a5bbbfc6871d044b2a09a78e963"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@cache
def prior_checks():
    if sha(prior.REPORT) != PRIOR_SHA or sha(remainder.REPORT) != A7_SHA:
        raise ValueError("A pinned A.9 or A.7 response input changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {"A9": PRIOR_SHA, "A7_calibration_input": A7_SHA}


def serialize(value):
    if isinstance(value, dict):
        return {name: serialize(item) for name, item in value.items()}
    if isinstance(value, list):
        return [serialize(item) for item in value]
    if isinstance(value, bool):
        return value
    return str(value)


def exact_checks():
    groups = {"actual_Dyson_and_Abel_response": mode_lipschitz.identities(),
              "causal_log_pole_cut_and_polynomial_data": log_inverse.identities()}
    for name, values in groups.items():
        if any(sp.simplify(value) != 0 for value in values.values()):
            raise ValueError("A response/inverse identity failed: "+name)
    return {name: dict.fromkeys(values, "0") for name, values in groups.items()}


def checked_constants():
    actual_modes, actual_inverse = mode_lipschitz.calibration(), log_inverse.calibration()
    independent_report = independent.replay(json.loads(remainder.REPORT.read_text()))
    mode_bridge = {
        "prepared_U1": actual_modes["prepared_U1"],
        "derivative_cap": actual_modes["full"]["derivative_cap"],
        "strength": actual_modes["full"]["strength"],
        "full_derivative": actual_modes["full"]["derivative"],
        "full_value": actual_modes["full"]["value"],
        "log_ratio_upper": actual_modes["shared"]["log_ratio_upper"],
        "shared_quadratic": actual_modes["shared"]["quadratic_derivative"],
        "shared_higher": actual_modes["shared"]["higher_derivative"],
        "shared_derivative": actual_modes["shared"]["derivative"],
        "shared_value": actual_modes["shared"]["value"],
        "full_rounding_gap": actual_modes["full_rounding_gap"],
        "shared_rounding_gap": actual_modes["shared_rounding_gap"],
    }
    if serialize(mode_bridge) != independent_report["mode_calibration"]:
        raise ValueError("The independent actual-response calibration disagrees")
    if any(actual_modes[key] <= 0 for key in ("full_rounding_gap", "shared_rounding_gap")):
        raise ValueError("A claimed response rounding gap is not positive")
    inverse = actual_inverse["isolated_inverse"]
    abstract = actual_inverse["abstract_conditional_example"]
    inverse_bridge = {key: inverse[key] for key in (
        "length", "split", "pole_norm_upper", "low_cut_norm_upper",
        "high_cut_norm_upper", "norm_upper")}
    inverse_bridge.update({"rounded_gap": sp.Rational(1, 3)-inverse["norm_upper"],
                           "abstract_q": abstract["q_upper"],
                           "abstract_distance": abstract["distance_upper"],
                           "abstract_ball_margin": abstract["self_map_margin"]})
    if serialize(inverse_bridge) != independent_report["inverse_calibration"]:
        raise ValueError("The independent causal-inverse calibration disagrees")
    for order in range(3, 10):
        symbolic = mode_lipschitz.dyson_coefficients(order, 3, sp.Rational(1, 3))
        bridge = {key: symbolic[key] for key in ("infrared", "ultraviolet", "total")}
        if serialize(bridge) != independent_report["dyson_allocation_checks"][str(order)]:
            raise ValueError("An independent two-mode Dyson allocation disagrees")
    return {"actual_mode_ball": serialize(actual_modes),
            "isolated_inverse_and_conditional_example": serialize(actual_inverse),
            "independent_Fraction_replay": independent_report}


def build_report():
    inputs = prior_checks()
    identities, constants = exact_checks(), checked_constants()
    paths = sorted(ROOT.glob("src/p8a_existence/*.py"))+sorted(ROOT.glob("tests/*.py"))
    paths += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1, "claim": "P8-A.10", "date": "2026-09-06",
        "status": "ACTUAL_NONLINEAR_C1_RESPONSE_AND_ISOLATED_CAUSAL_INVERSE_CERTIFIED; FULL_SEE_EXISTENCE_AND_P8_OPEN",
        "prior_sha256": inputs,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in paths},
        "formulation": "FORMULATION.md",
        "written_proofs": ["notes/mode-response.md", "notes/inverse.md"],
        "source_audit": "notes/sources.md", "remaining_obligations": "notes/see-gap.md",
        "exact_residuals": identities, "derived_constants": constants,
        "mode_theorem": {
            "domain": "real U,V in C1[0,T], U(0)=V(0)=0, common original plane-wave data; no sign or mass-gap requirement",
            "bounds": "||U'||,||V'||<=M; D=||U'-V'||, z=M*T^3",
            "functional": "R=integral k*(|w_k|^2-1-F1_k) dk; quadratic Abel limit plus absolutely convergent higher Dyson orders",
            "actual_state_scope": "agrees with the actual transported Hadamard state's nonlinear Wick term for smooth past-flat geometries; generic C1 inputs are not asserted Hadamard",
            "full_derivative_coefficient": "2*z+18*z^2*exp(2*z)",
            "full_value_coefficient": "T*(2*z+18*z^2*exp(2*z))",
            "shared_history": "U=V on [0,T-L], 0<L<=T; retain original nonzero prepared history",
            "shared_derivative_coefficient": "M*T*L^2*(5/4+log(T/L)/2)+18*M^2*T^5*L*exp(2*M*T^3)",
            "shared_value_coefficient": "L times shared derivative coefficient",
            "higher_potential_derivatives_needed": False,
            "analytic_strength_domain": "all finite M>=0,T>0",
            "rational_API_strength_domain": "0<=M*T^3<=1/4",
            "actual_calibration": "dimensionless delta=10^-14, T<=3, M=2*delta*b1 with b1 from pinned A7, future L=1/4; derivative coefficients <10^-8 and <10^-10",
            "calibration_clock": "eta_star is the last free slice; the original state is freely propagated from y=1/2; T is only a history upper bound, not an unproved future extension",
            "full_SEE_Lipschitz_constant": False,
            "uniqueness_from_density_of_smooth_past_flat_subset": False,
        },
        "inverse_theorem": {
            "operator": "D_beta f=[(beta-EulerGamma)*(f-f(0))-integral_0^t log(t-s)*f'(s) ds]/(8*pi^2)",
            "Laplace_multiplier": "(beta+log(s))/(8*pi^2) acting on f-f(0)",
            "causal_kernel": "8*pi^2*[p*exp(p*t)+integral_0^infinity exp(-r*t)/((log(r)+beta)^2+pi^2) dr], p=exp(-beta)>0",
            "exact_C0_norm": "8*pi^2*[exp(p*L)-1+integral_0^infinity (1-exp(-r*L))/(r*((log(r)+beta)^2+pi^2)) dr]",
            "finite_upper_bound": "8*pi^2*[exp(p*L)-1+sqrt(L)/pi^2+1/(beta-log(L)/2)] for beta-log(L)/2>0",
            "pole_omitted": False,
            "C1_endomorphism_or_automatic_C1_inverse": False,
            "abstract_dyadic_example": "beta>=-1, L=2^-1024 gives norm<1/3; not a physical proper-time or SEE interval",
            "finite_stress_prescription_changed": False,
            "full_SEE_map_identified_with_isolated_operator": False,
        },
        "conditional_fixed_point_gate": {
            "necessary_inputs": "proved concrete map on a closed positive-radius C0 ball; full Lipschitz Lambda, operator bound C, fixed-point residual eta, radius r",
            "verified_inequalities": "q=C*Lambda<1 and eta+q*r<=r; unique ball fixed point then has error <=eta/(1-q)",
            "illustrative_inputs_are_actual_SEE_data": False,
            "domain_or_input_bounds_proved_by_numeric_API": False,
        },
        "negative_controls": {
            "nonzero_initial_potential": "derivative transfer acquires -sin(k*t)*U(0)/k at first order",
            "quadratic_absolute_UV_majorant": "integral k^-1 diverges; exact two-mode oscillatory/Abel cancellation is required",
            "discarding_common_prehistory": "would lose the future-past mixed quadratic and marked Dyson terms",
            "missing_pole_Laplace_defect": str(log_inverse.controls()["omitted_pole_Laplace_defect"]),
            "C1_domain_control": "D_beta[t] has an unbounded initial derivative and J_beta[1] is not C1 at zero",
            "inexact_or_outside_domain_inputs": "rejected, including booleans, binary floats, nonfinite values, forbidden order, nonpositive duration/split, and failed strict contraction",
        },
        "verification_boundary": [
            "Pinned A9 and its A8/A7 lineage are replayed without old-file edits",
            "New constants have independently assembled Fraction-only checks; no sampled proof of an infinite-frequency claim",
            "Independent actual quadratic product, radial measure, Abel, initial endpoint and clock-normalization tests are hashed",
            "Integral convergence, keyhole contour, marked simplex and Banach implications are explicit written proofs, not formalized PDE results",
            "Independent review corrected the C1 extension-density wording before the source hashes were frozen",
        ],
        "not_established": [
            "Full actual SEE functional self-map, gravitational/geometric Lipschitz constants or long-interval shadowing",
            "Exact compatible initial Einstein constraint with fixed radiation and a joint metric/state preparation",
            "A smooth Hadamard exact solution, or transfer of A9's fixed-metric QSEI to such a solution",
            "An arbitrary-state, massive/interacting/realistic-field response theorem",
            "A cosmological incompleteness theorem or completion of P8a/P8",
        ],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Actual-response/inverse certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8(a) A.10: actual C1 mode response and causal inverse replay passed; full SEE and P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
