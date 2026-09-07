"""Read-only determinant-sum background theorem and source-pinned replay."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_star import verify as prior
from p8_trimetric import verify as action_prior

from . import background, controls, independent, potential

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"determinant-no-bounce.json"
P8 = next(path for path in ROOT.parents if path.name == "P8")
PRIOR_SHA = "54f23dfd05d7580969550d06ac4a72adc5faa61988be9e7c95eeb199fb546353"
ACTION_SHA = "062ea71b4fcb139af3ea727f013eeeb3a340b6e37207cb872e372eacce09631e"
CONTRACT_SHA = "d00df35d5821b05baa61e897725a22ed1ac21d0536b74bdaece3ab0093215901"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@cache
def prior_checks():
    pins = {"S6_16_frozen_STAR_and_arithmetic_lineage": (prior.REPORT, PRIOR_SHA),
            "S6_13_literal_auxiliary_action": (action_prior.REPORT, ACTION_SHA),
            "adopted_S6_matching_contract": (P8/"s6"/"FORMULATION.md", CONTRACT_SHA)}
    if any(sha(path) != expected for path, expected in pins.values()):
        raise ValueError("A frozen action, arithmetic or matching-contract input changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {key: expected for key, (_, expected) in pins.items()}


def bridges():
    stresses, rates = [], []
    for fixture in independent.stress_fixtures():
        actual = potential.evaluate(fixture["beta"], fixture["a"], fixture["n"], lambda_=fixture["lambda"])
        keys = ("S_a", "S_N", "rho", "pressure", "nulls")
        if any(actual[key] != fixture[key] for key in keys):
            raise ValueError("Primary versus full-coframe Fraction variation failed")
        stresses.append({"regular_sum": actual["regular_sum"], "active": list(actual["active"]),
                         "primary_minus_independent": dict.fromkeys(keys, "0")})
    for fixture in independent.rate_fixtures():
        actual = background.reconstruct(fixture["Gs"], fixture["ys"], fixture["cs"],
                                        fixture["H"], fixture["nulls"])
        keys = ("K", "Kprime", "Hprime", "yprimes", "weighted_null")
        if any(actual[key] != fixture[key] for key in keys):
            raise ValueError("Primary versus independent physical-clock/source weighting failed")
        if actual["scaled_H_prime"].is_nonpositive is not True:
            raise ValueError("The positive-source dynamic sign was lost")
        rates.append({"active_count": len(fixture["Gs"]), "primary_minus_independent": dict.fromkeys(keys, "0")})
    return {"full_coframe_variation": stresses, "proper_clock_and_null_sources": rates}


def omission_controls():
    singular = controls.singular_sum_bianchi_control()
    if (singular["S_a"], singular["S_N"], singular["S_a_dot"]) != (2, 0, 0):
        raise ValueError("The singular-sum Bianchi-only countercontrol changed")
    if any(singular["bianchi"]) or len(set(singular["Q_i"])) != 3:
        raise ValueError("Unjustified singular-sum velocity locking was not excluded")
    if singular["actual_solution"] or singular["branch_health_claim"]:
        raise ValueError("A Bianchi-only point was promoted to an actual healthy solution")
    broken = potential.evaluate([1, -1], [1, 1], [1, 2])
    if broken["regular_sum"] or broken["rho"] != (0, 0) or broken["pressure"] != (0, 0):
        raise ValueError("The spatial-sum-zero polynomial control failed")
    moved_source = controls.source_relocation_control()
    if moved_source != (-sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(1, 2)):
        raise ValueError("Moving matter to the auxiliary field was silently allowed")
    weighted = background.reconstruct([1, 2], [1, 3], [1, 2], 0, [1, 2])
    if weighted["K"] != 19 or weighted["weighted_null"] != 325:
        raise ValueError("A physical inverse-ratio or source-measure omission control failed")
    rejected = 0
    bad_calls = [lambda: potential.evaluate([1.0], [1], [1]),
                 lambda: potential.evaluate([1], [True], [1]),
                 lambda: potential.evaluate([1], [1], [0]),
                 lambda: potential.evaluate([1], [1], [1], lambda_=0),
                 lambda: background.reconstruct([sp.oo], [1], [1], 0, [0]),
                 lambda: background.reconstruct([1], [1], [1], 0, [-1]),
                 lambda: background.reconstruct([1], [2], [1], 0, [0]),
                 lambda: background.reconstruct([], [], [], 0, [])]
    for call in bad_calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(bad_calls):
        raise ValueError("An invalid theorem-domain numeric input was admitted")
    return {"singular_time_sum": {"S_a": "2", "S_N": "0", "dot_S_a": "0",
                                  "unlocked_velocities": ["1", "-1", "0"], "actual_solution": False},
            "spatial_zero_sum_interaction_stresses": "all zero; not a classification of S_N-zero branches",
            "canonical_source_moved_to_w_Euler": [str(value) for value in moved_source],
            "correct_K": "19", "wrong_inverse_ratio_K": "11/9",
            "correct_weighted_source": "325", "wrong_unweighted_source": "3",
            "invalid_domain_calls_rejected": rejected}


@cache
def build_report():
    previous = prior_checks()
    residuals = {"literal_determinant_and_unfactored_Bianchi": potential.checks(),
                 "actual_clock_weighted_Einstein_identity": background.checks(),
                 "actual_solutions_auxiliary_map_and_countercontrols": controls.checks()}
    if any(sp.simplify(value) != 0 for group in residuals.values() for value in group.values()):
        raise ValueError("A determinant-action or actual-background identity failed")
    files = sorted(ROOT.glob("src/p8_determinant/*.py"))+sorted(ROOT.glob("tests/*.py"))
    files += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1, "claim": "P8-S6.18.DETERMINANT", "date": "2026-09-06",
        "status": "REGULAR_ALIGNED_FLAT_DETERMINANT_SUM_SEPARATE_NEC_SOURCES_NO_BOUNCE; MATCHING_AND_P8_B_OPEN",
        "prior_context_sha256": previous,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in files},
        "formulation": "FORMULATION.md", "written_proof": "notes/proof.md", "source_audit": "notes/sources.md",
        "exact_residuals": {name: dict.fromkeys(group, "0") for name, group in residuals.items()},
        "independent_Fraction_replay": independent.checks(),
        "primary_independent_bridges": bridges(), "checked_omission_controls": omission_controls(),
        "fixed_action_domain": {
            "action": "-sum_i G_i/2 sqrt|g_i|(R_B+2Lambda_i)-lambda det(sum_i beta_i e_i)+sum_i S_m,i[g_i,psi_i]",
            "constants": "lambda=m^4>0; finite G_i>0; real constant beta_i,Lambda_i; finite nonempty metric set",
            "geometry": "Common aligned diagonal positive smooth flat-FLRW coframes, connected open chart; not every relative-Lorentz branch",
            "regular_sum": "S_a=sum beta_i a_i !=0 and S_N=sum beta_i n_i !=0 everywhere in the theorem interval",
            "signs": "Same-sign nonzero beta is a sufficient source-paper specialization, not necessary in the algebraic proof; mixed signs get no health verdict",
            "physical_sources": "Independent on-shell homogeneous/isotropic sources each minimally coupled to its own actual EH metric, separately conserved and NEC; no matter on summed or auxiliary metric",
            "curvature_and_clock": "+---,R_B=-6(DH+2H²); 3G_i H_i²=G_i Lambda_i+rho_m,i+rho_int,i; physical dT=n_r dt",
        },
        "theorem": {
            "interaction": "rho_int,i=lambda beta_i S_a³/a_i³; p_int,i=-lambda beta_i S_N S_a²/(n_i a_i²)",
            "literal_null_sum": "sum_i n_i a_i³ (rho_int,i+p_int,i)=0, without dividing by either sum",
            "unfactored_Bianchi": "3lambda beta_i S_a²/(n_i a_i³)[dot S_a-S_N dot a_i/n_i]=0",
            "locking_domain": "For active beta_i!=0 and regular sum, dot a_i/n_i=dot S_a/S_N even at simultaneous zero velocities",
            "kinematics": "y_i=a_i/a_r,c_i=n_i/(n_r y_i),y_r=c_r=1; H_i=H/y_i,y_i'=y_i(c_i-1)H",
            "positive_weight": "K=sum_active G_i y_i²>0; K'=2H sum_active G_i y_i²(c_i-1)",
            "exact_identity": "K H'-H K'/2=-sum_active c_i y_i^4(rho_m,i+p_m,i)/2",
            "monotonicity": "(H/sqrtK)'=-sum_active c_i y_i^4(rho_m,i+p_m,i)/(2K^(3/2))<=0",
            "conclusion": "H(T0)<=0 implies H(T)<=0 for every later T in the same regular interval; no degenerate or ordinary contraction-to-expansion",
            "zero_beta": "Omit disconnected fields from active K; if the physical beta_r=0, its separate flat GR H'=-nu_r/(2G_r)<=0 applies without assumptions on unrelated sums",
        },
        "actual_controls": {
            "rolling": "lambda=1,beta=(1,-1,1),G=(1,2,3),y=(1,2,3),Lambda=(-8,1/2,-8/81); a_i=y_i T^(1/3),n_i=y_i,phi_i=sqrt(2G_i/3)log(T/T*)",
            "rolling_Einstein_and_scalar": "Each rho_m,i=p_m,i=G_i/(3y_i²T²); S_N=2,S_a=2T^(1/3),K=36,normalized derivative=-1/(18T²); all individual Einstein and scalar equations vanish",
            "vacuum": "Same beta,G,y; a_i=y_i exp(H0 T),n_i=y_i,zero matter,Lambda_i=3H0²/y_i²-rho_int,i/G_i; normalized derivative zero for any fixed real H0",
            "meaning": "Actual regular backgrounds, not a bounce, perturbative-health or original C/D matching witness",
        },
        "auxiliary_source_preserving_identity": {
            "map": "L_aux=lambda det w(3-tr(w^-1 U)),B=-3lambda/2,p_i=lambda beta_i/2; invertible stationary w=U gives -lambda detU",
            "source_location": "Matter remains on original EH leaves; adding matter on w changes its Euler equation and generally its stationary value",
            "chart": "A positive auxiliary-w chart only where applicable; the direct determinant proof needs actual positive EH coframes and nonzero sums, not positive w",
        },
        "verification_boundary": "Literal symbolic, separate Fraction/full-coframe/adjugate replay and reserved independent audit support the written monotonicity proof; not proof-assistant formalized",
        "primary_sources": ["https://arxiv.org/pdf/2607.24347", "https://arxiv.org/pdf/1804.04671"],
        "not_established": [
            "A theorem through persistent S_N=0,S_a!=0 or any classified singular-sum continuation",
            "Actual solution or perturbative health of the Bianchi-only singular countercontrol",
            "Uniqueness of the common-Lorentz branch or a no-bounce theorem for nonaligned coframes",
            "An unchanged theorem after moving matter to the composite/auxiliary metric, using variable/derivative interactions, nonflat slices or physical NEC violation",
            "Scalar/vector/tensor health, a quantum cutoff, Regge/loop control or a full/light-EFT UV verdict from the paper's high-density cone speeds",
            "Original C/D matching, a universal parent exclusion, or closure of S6, P8(b), or P8",
        ],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Determinant-sum no-bounce certificate differs from exact read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.18.DETERMINANT: regular-sum physical no-bounce replay passed; singular sums and P8(b) OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
