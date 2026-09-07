"""Read-only robust cosmological-strength photon theorem calibration replay."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8a_maxwell_focusing import verify as prior

from . import calibration, controls, dictionary, geometry, independent

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"cosmological-calibration.json"
PARENT_SHA = "22cfba547368b420767bbe87068c5b45fb556992c9e37a389091525f3e2ef62a"


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
        raise ValueError("the pinned photon history-focusing theorem changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {"A17": PARENT_SHA}


def exact_checks():
    groups = {"all_p_proper_clock_and_history": geometry.identities(),
              "finite_beta_and_sourced_focusing_margin": calibration.identities(),
              "physical_clock_and_signed_source_dictionary": dictionary.identities(),
              "complete_geometric_family_and_tight_envelope_exclusion": controls.identities()}
    for name, group in groups.items():
        if any(sp.simplify(value) != 0 for value in group.values()):
            raise ValueError("an exact cosmological-strength identity failed: "+name)
    return {name: dict.fromkeys(group, "0") for name, group in groups.items()}


def checked_constants():
    # A17's full replay above checks A16's exact report and source lineage.
    # The transitive photon report supplies monomials, not a second unpinned
    # field/prescription input.
    result = independent.replay(json.loads(prior.REPORT.read_text()),
                                json.loads(prior.prior.REPORT.read_text()))
    groups = {"geometry": geometry.calibration(), "calibration": calibration.calibration(),
              "dictionary": dictionary.calibration_data(), "controls": controls.calibration()}
    for name, group in groups.items():
        if serialize(group) != result[name]:
            raise ValueError("independent Fraction cosmology replay disagrees: "+name)
    return {**{name: serialize(group) for name, group in groups.items()},
            "independent_Fraction_replay": result}


def build_report():
    inputs = prior_checks()
    identities, constants = exact_checks(), checked_constants()
    paths = sorted(ROOT.glob("src/p8a_maxwell_cosmology/*.py"))+sorted(ROOT.glob("tests/*.py"))
    paths += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1, "claim": "P8-A.18", "date": "2026-09-06",
        "status": "ROBUST_FREE_PHOTON_FLAT_FLRW_COSMOLOGICAL_STRENGTH_INCOMPLETENESS_CALIBRATION_CERTIFIED; ROOT_SCOPE_ASSESSMENT_SEPARATE",
        "prior_sha256": inputs,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in paths},
        "formulation": "FORMULATION.md", "written_proof": "notes/proof.md",
        "source_audit": "notes/sources.md", "literal_objective_assessment": "notes/closure.md",
        "exact_residuals": identities, "derived_constants": constants,
        "theorem": {
            "spacetime": "stipulated smooth positive global spatially flat FLRW product I_s times R^3",
            "direction": "s is the chosen contracting proper normal parameter; past incompleteness follows with s=T_star-T in an expanding interpretation",
            "actual_field": "one physical free Maxwell field, arbitrary Hadamard state, no added state symmetry or quasifree restriction",
            "actual_source_equation": "G_FK+Lambda g_FK=-kappa(T_Maxwell+T_other), E_other>=ell",
            "prescription": "same A16 conserved real beta_M R^2 family, separately fixed Lambda/Newton coupling; no inherited scalar gamma",
            "actual_history": "[-tau/100,0] including smooth endpoint neighborhoods actually exists and lies in the stated C3 tube",
            "reference_family": "H_p(x)=-p/(p/2-x), p in [1/2,2/3], x=s/tau in [-1/100,0]",
            "C3_error_caps": ["1/100", "1/2", "3", "16"],
            "actual_past_contraction": "H<=-(19/10)/tau",
            "actual_past_caps": ["21/10", "9", "70", "800"],
            "future_caps_conditional_on_reaching_tau": ["4", "128", "16384", "1048576"],
            "future_cap_form": "|H^(j)(s)|<=c_j/(tau-s)^(j+1),0<=s<tau,j0..3, only if a normal reaches tau",
            "future_caps_derived_from_history_or_QSEI": False,
            "generic_quantum_gate": "delta(1+abs(beta_M))<=1e-8, delta=kappa*hbar/(8*pi^2*tau^2)",
            "source_gate": "sigma=tau^2[max(Lambda,0)+kappa max(-ell,0)]<=5",
            "strict_margin_lower": "47079/350000 >1/8",
            "conclusion": "upper endpoint of I_s<=tau; all directed timelike curve lengths from s0<=tau; stipulated spacetime timelike geodesically incomplete",
            "initial_SEC_is_an_independent_energy_premise": False,
            "this_small_calibration_tube_has_initial_timelike_convergence": True,
            "future_pointwise_SEC_imposed": False,
        },
        "observer_and_physical_dictionary": {
            "reference_Hstar_times_tau": "2",
            "actual_H0_times_tau_unanchored": ["199/100", "201/100"],
            "optional_anchored_slice": "E(0)=0 makes actual H0*tau=2; relative fixed-endpoint neighborhood only",
            "scale_factor_comparison": "normalized a(s)/a(0) versus normalized reference, or matched a(0)",
            "reference_age_not_actual_quantum_age": "T_star/tau=p/2 in[1/4,1/3]",
            "seconds_Planck_gate": "tau_seconds/t_P>=1e4*sqrt((1+abs(beta_M))/pi)",
            "reference_Lambda_fraction_upper": "5/12",
            "unanchored_actual_Lambda_fraction_upper": "50000/118803",
            "pure_Maxwell_specialization": "T_other=0,ell=0; Lambda remains separately specified",
            "ordinary_matter_corollary": "nonnegative radiation/dust EED requires explicit separate rho>=0 and equation-of-state assumptions",
            "observed_parameters_or_actual_quantum_comparator_used": False,
        },
        "geometric_nonvacuity_controls": {
            "complete_family": "all p in[1/2,2/3]: compact future P7 cutoff at1/8, positive unit-L1 mollification width1e-12",
            "raw_regularity": "C3 and W4,infinity before mollification; smooth afterwards",
            "future_timelike_and_null_complete": True,
            "satisfies_actual_past_tube_and_final_future_caps": True,
            "actual_allowed_Maxwell_SEE_solution": False,
            "geometric_assumptions_alone_force_incompleteness": False,
            "rejected_tight_future": "c=d gives a Taylor contradiction atx1/20 with gap4707907/62554080 before any QSEI",
        },
        "verification_boundary": [
            "A17 report and full source-hashed lineage replayed read-only before importing the A16 photon coefficient polynomial",
            "Continuous all-p/all-history bounds and smooth noncompact complete comparator proved analytically, not sampled numerically",
            "Independent stdlib Fraction polynomial and pinned-coefficient reconstruction checks the complete core calibration",
            "Parent-owned covariant/time-orientation/Leibniz/source and omission audit included unchanged",
            "Machine algebra and source hashes support written proofs; no proof-assistant formalization is claimed",
        ],
        "scope_assessment": {
            "proved_specialization": "free-photon global spatially flat FLRW incompleteness with realistic-field and cosmological-strength constants under explicit geometry and source assumptions",
            "actual_SEE_witness_required_for_conditional_theorem": False,
            "unrestricted_field_or_spacetime_theorem": False,
            "root_P8_completion_declared": False,
        },
        "not_established": [
            "observed-universe parameter fit or an actual quantum SEE solution equal to a classical power-law/comparison metric",
            "conditional future curvature bounds derived from the QSEI or the actual short history",
            "a quantum energy inequality for added ordinary matter without its own explicit assumptions",
            "interacting QED, arbitrary realistic fields, null/boosted focusing or arbitrary inhomogeneous spacetimes",
            "curvature blow-up, maximal inextendibility, optimal constants or fundamental EFT validity at an endpoint",
            "root-level closure of P8a, P8b or P8 as a whole",
        ],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("the robust photon cosmology certificate differs from exact read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8(a) A.18: robust photon flat-FLRW cosmological-strength incompleteness calibration passed; root scope assessment separate")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
