"""Read-only actual thermal-photon state, SEE history and endpoint replay."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8a_maxwell_cosmology import verify as prior

from . import bounds, dynamics, independent, state

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"thermal-see.json"
PARENT_SHA = "b03767fb43c49a0f6c5ff34307355c4cd37a75930f61e07f4826c06101e4d23d"
PHOTON_SHA = "f59d63524bf1d98998a928fbe5139eef45f3703ebd1e61cbdf470e738ce94b94"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serialize(value):
    if isinstance(value, dict):
        return {key: serialize(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(item) for item in value]
    if isinstance(value, bool):
        return value
    return str(value)


@cache
def prior_checks():
    if sha(prior.REPORT) != PARENT_SHA:
        raise ValueError("the pinned A18 cosmological calibration changed")
    if sha(prior.prior.prior.REPORT) != PHOTON_SHA:
        raise ValueError("the transitive pinned A16 photon input changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {"A18": PARENT_SHA, "A16_transitive": PHOTON_SHA}


def exact_checks():
    groups = {"physical_thermal_state_and_prescription": state.identities(),
              "exact_SEE_clock_and_low_branch": dynamics.identities(),
              "continuous_C3_response_and_domain": bounds.identities()}
    for name, group in groups.items():
        if any(sp.simplify(value) != 0 for value in group.values()):
            raise ValueError("an exact thermal-photon identity failed: "+name)
    return {name: dict.fromkeys(group, "0") for name, group in groups.items()}


def checked_constants():
    # prior_checks first replays A18 and its complete A17/A16 source lineage.
    # No unpinned alternative photon prescription is used by this engine.
    prior_checks()
    result = independent.replay(json.loads(prior.REPORT.read_text()),
                                json.loads(prior.prior.prior.REPORT.read_text()))
    groups = {"state": state.calibration(), "dynamics": dynamics.calibration(),
              "bounds": bounds.calibration()}
    for name, group in groups.items():
        if serialize(group) != result[name]:
            raise ValueError("independent Fraction thermal replay disagrees: "+name)
    return {**{name: serialize(group) for name, group in groups.items()},
            "independent_Fraction_replay": result}


def build_report():
    inputs = prior_checks()
    identities, constants = exact_checks(), checked_constants()
    paths = sorted(ROOT.glob("src/p8a_maxwell_thermal/*.py"))+sorted(ROOT.glob("tests/*.py"))
    paths += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1, "claim": "P8-A.19", "date": "2026-09-06",
        "status": "ACTUAL_THERMAL_PHOTON_NAMED_SEE_HISTORY_AND_FINITE_LOW_BRANCH_ENDPOINT_CERTIFIED; OPTIONAL_STRENGTHENING",
        "prior_sha256": inputs,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in paths},
        "formulation": "FORMULATION.md", "written_proof": "notes/proof.md",
        "source_audit": "notes/sources.md", "exact_residuals": identities,
        "derived_constants": constants,
        "actual_state": {
            "field": "one physical free Maxwell field with two transverse polarizations",
            "construction": "positive quasifree occupation state, n(k)=(exp(b_T*k)-1)^-1, then gauge-invariant conformal transport",
            "Hadamard": "thermal-minus-vacuum field-strength covariance is smooth by UV exponential and IR integrable bounds",
            "domain": "every interior point of the selected smooth positive flat-FLRW open strip",
            "thermal_amplitude": "Q=hbar*pi^2/(15*b_T^4)",
            "physical_temperature": "k_B*T_physical=hbar/(a*b_T), c=1",
            "b_T_is_photon_beta_M_or_scalar_gamma": False,
            "global_cosmic_time_equilibrium_claimed": False,
            "every_Hadamard_state_solves_the_same_metric": False,
        },
        "named_source_and_solution": {
            "prescription": "A16 beta_M=0, Lambda=0, no additional matter or independent nonzero gravitational curvature-squared coupling",
            "equation": "G_FK=-kappa*T_Maxwell, kappa>0,hbar>0",
            "physical_density": "Q/a^4+31*hbar*H^4/(480*pi^2)",
            "exact_constraint": "H^2-b*H^4=kappa*Q/(3*a^4), b=31*kappa*hbar/(1440*pi^2)",
            "independent_pressure_equation_checked": True,
            "additional_classical_radiation_constant_assumed": False,
            "vacuum_H4_anomaly_omitted": False,
            "order_reduced_surrogate_or_only_formal_solution": False,
            "delta": "kappa*hbar/(8*pi^2*tau^2), 0<delta<=1e-8",
            "lambda": "b/tau^2=31*delta/180",
            "normal_clock": "s=T0-T, x=s/tau, y=-tau*H_normal; expansion-oriented interpretation has a past endpoint",
            "anchor": "a(0)=1,y(0)=2,Q=12*(1-4*lambda)/(kappa*tau^2)",
            "state_temperature_chosen_from_constraint": "b_T=(hbar*pi^2/(15*Q))^(1/4)",
            "admitted_branch": "0<y<1/sqrt(2*lambda); positive algebraic a^4 on another branch is insufficient",
            "exact_clock": "x=Phi(2)-Phi(y), Phi(y)=1/(2*y)+sqrt(lambda)*atanh(sqrt(lambda)*y)/2",
            "exact_scale": "a^4=4*(1-4*lambda)/(y^2*(1-lambda*y^2))",
        },
        "actual_history_embedding": {
            "interval": "all x in[-1/100,0], with smooth endpoint neighborhoods",
            "reference": "anchored A18 radiation p=1/2, y_rad=2/(1-4*x)",
            "continuous_C3_coefficients_per_lambda": ["16/25", "1728/25", "155136/25", "14770176/25"],
            "strictly_inside_pinned_error_caps": ["1/100", "1/2", "3", "16"],
            "actual_observer_H0_times_tau": "2",
            "absolute_scale_normalization_inferred_from_jet_closeness": False,
            "scale_a0_one_is_explicit_choice": True,
            "actual_SEE_history_and_state_realization": True,
            "all_radiation_to_matter_p_values_realized": False,
            "continuous_proof_not_a_sample_scan_or_unknown_asymptotic_remainder": True,
            "sigma_additional_source_budget": "0",
        },
        "endpoint": {
            "regular_open_x_domain": "(-infinity,x_end)",
            "exact_x_end": "Phi(2)-Phi(1/sqrt(2*lambda))",
            "strict_x_end_bounds": ["1/16", "1/4"],
            "a_end_fourth": "16*lambda*(1-4*lambda)>0 for each fixed delta>0",
            "uniform_positive_a_end_as_delta_tends_to_zero": False,
            "Ricci_scalar": "tau^2*R_FK=-12*lambda*y^4/(1-2*lambda*y^2) tends to -infinity",
            "endpoint_admitted_as_smooth_metric_or_state_initial_data": False,
            "conclusion": "selected actual low-branch spacetime is timelike geodesically incomplete; finite comoving remaining proper length tau*x_end",
            "C2_continuation_through_same_endpoint_metric": False,
            "fundamental_EFT_validity_at_anomaly_scale_proved": False,
            "A18_future_samplers_or_caps_evaluated_through_endpoint": False,
        },
        "relationship_to_QSEI_objective": {
            "particular_branch_has_positive_EED": True,
            "new_QEI_only_endpoint_proof": False,
            "A18_all_Hadamard_theorem_changed": False,
            "actual_history_realization_is_optional_strengthening": True,
            "reopens_scoped_P8a_completion": False,
            "root_P8_status_assigned_by_this_child": False,
        },
        "verification_boundary": [
            "A18 and its full immutable A17/A16 lineage replayed read-only before photon coefficients are imported",
            "Actual state positivity/Hadamard property and global open-branch continuation proved analytically, not certified by numerical sampling",
            "Exact symbolic source/clock/jet identities supplemented by independent Fraction reconstruction and parent-owned covariant audit",
            "Proof and arithmetic certificate are not proof-assistant formalization or external peer review",
        ],
        "not_established": [
            "actual thermal SEE existence for generic beta_M, nonzero Lambda, added matter or independent curvature-squared couplings",
            "all-Hadamard target states satisfying the same selected thermal background equation",
            "interacting QED or observationally fitted cosmic histories",
            "fundamental EFT control at the finite high-curvature endpoint or smooth passage to another algebraic branch",
            "a new QEI-only incompleteness proof or a new completion prerequisite for A18",
        ],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("the thermal-photon certificate differs from exact read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8(a) A.19: actual thermal-photon state, exact SEE radiation history and finite low-branch endpoint passed; optional strengthening")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
