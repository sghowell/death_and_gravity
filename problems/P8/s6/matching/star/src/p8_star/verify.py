"""Read-only exact constant-star no-bounce certificate and immutable replay."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_trimetric import verify as action_prior
from p8_trimetric_cones import verify as prior

from . import background, branches, independent, potential

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"star-no-bounce.json"
P8 = next(path for path in ROOT.parents if path.name == "P8")
PRIOR_SHA = "ab280acf8d20fb2ae18d75c265deecd59832e993ba0338e9e5c8f7691c0a7e6e"
ACTION_SHA = "062ea71b4fcb139af3ea727f013eeeb3a340b6e37207cb872e372eacce09631e"
CONTRACT_SHA = "d00df35d5821b05baa61e897725a22ed1ac21d0536b74bdaece3ab0093215901"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@cache
def prior_checks():
    pins = {"S6_14_physical_u_and_curvature_lineage": (prior.REPORT, PRIOR_SHA),
            "S6_13_literal_auxiliary_action": (action_prior.REPORT, ACTION_SHA),
            "adopted_S6_conditional_matching_contract": (P8/"s6"/"FORMULATION.md", CONTRACT_SHA)}
    if any(sha(path) != expected for path, expected in pins.values()):
        raise ValueError("A pinned action, convention or matching contract changed")
    # S6.14 recursively validates S6.13 and its arithmetic/action ancestry.
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {key: expected for key, (_, expected) in pins.items()}


def independent_bridges():
    stress = []
    for fixture in independent.stress_fixtures():
        actual = potential.evaluate(fixture["beta"], fixture["R"], fixture["N"])
        keys = ("rho_i", "pressure_i", "rho_u", "pressure_u", "null_i", "null_u")
        if any(actual[key] != fixture[key] for key in keys):
            raise ValueError("Literal independent coframe jet versus primary variation failed")
        stress.append({"R": str(fixture["R"]), "N": str(fixture["N"]),
                       "primary_minus_independent": dict.fromkeys(keys, "0")})
    dynamics = []
    for fixture in independent.dynamic_fixtures():
        actual = background.dynamic_reconstruction(
            fixture["Gs"], fixture["Rs"], fixture["cs"], fixture["H"], fixture["nh"],
            central_G=fixture["G_u"])
        keys = ("K", "Kprime", "Hprime", "Rprimes")
        if any(actual[key] != fixture[key] for key in keys):
            raise ValueError("Independent finite-star chain rule versus primary rate formula failed")
        if actual["scaled_H_prime"].is_nonpositive is not True:
            raise ValueError("The conditional dynamic-stratum sign failed")
        dynamics.append({"leaf_count": len(fixture["Gs"]), "G_u": str(fixture["G_u"]),
                         "primary_minus_independent": dict.fromkeys(keys, "0")})
    return {"literal_coframe_stresses": stress, "independent_clock_rate_fixtures": dynamics}


def controls():
    algebraic, disconnected = branches.algebraic_control(), branches.endpoint_control()
    t = algebraic["T"]
    if sp.simplify(algebraic["false_K_all_defect"]+1/(6*sp.sqrt(2)*t**2)) != 0:
        raise ValueError("The false all-leg monotonicity identity was not excluded")
    if algebraic["false_K_all_defect"].is_negative is not True:
        raise ValueError("The algebraic countercontrol must have a strict nonzero defect")
    if disconnected["H_u_prime_at_zero"] != 2 or any(
            disconnected[key] != 0 for key in ("central_Euler", "leaf_Euler", "scalar_Euler")):
        raise ValueError("The zero-central-kinetic disconnected exception failed")
    # All eight literal interaction variations and actual background equations
    # of the disconnected control are independently derived, not read as zeros.
    independent.disconnected_fixtures()
    beta = (1, 2, 3, 4, 5)
    correct_lapse = potential.evaluate(beta, 2, 3)["null_u"]
    wrong_lapse = potential.evaluate(beta, 2, sp.Rational(2, 3))["null_u"]
    with_center = background.dynamic_reconstruction([1], [1], [1], 0, 2, central_G=3)["Hprime"]
    without_center = background.dynamic_reconstruction([1], [1], [1], 0, 2)["Hprime"]
    if not correct_lapse > 0 > wrong_lapse or (with_center, without_center) != (-sp.Rational(1, 4), -1):
        raise ValueError("A proper-lapse or central-Einstein omission control failed")
    classifications = {
        "endpoint_only": branches.classify([7, 0, 0, 0, -5]),
        "old_nonzero_beta3": branches.classify([0, 0, 0, -2, 0]),
        "linear_positive_root": branches.classify([0, 1, "-1/2", 0, 0]),
        "double_positive_root": branches.classify([0, 1, -1, 1, 0]),
        "two_positive_roots": branches.classify([0, 2, "-3/2", 1, 0]),
        "irrational_positive_root": branches.classify([0, -2, 0, 1, 0]),
    }
    expected_roots = ((), (), ((1, 1),), ((1, 2),), ((1, 1), (2, 1)), ((sp.sqrt(2), 1),))
    if tuple(value["positive_roots"] for value in classifications.values()) != expected_roots:
        raise ValueError("A genuine-polynomial root classification control failed")
    rejected = 0
    calls = [lambda: potential.evaluate([0, 1.0, 0, 0, 0], 1, 1),
             lambda: potential.evaluate([0, 1, 0, 0, 0], True, 1),
             lambda: potential.evaluate([0, 1, 0, 0, 0], 1, 0),
             lambda: potential.evaluate([0, 1, 0, 0, 0], sp.oo, 1),
             lambda: background.dynamic_reconstruction([], [], [], 0, 0),
             lambda: background.dynamic_reconstruction([1], [1], [1], 0, -1),
             lambda: background.dynamic_reconstruction([1], [1], [1], 0, 1, central_G=-1)]
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An invalid theorem-domain input was silently admitted")
    return {"actual_algebraic_false_all_K_normalized_defect": str(algebraic["false_K_all_defect"]),
            "actual_disconnected_physical_Hprime_at_zero": "2",
            "correct_N_central_null": str(correct_lapse), "wrong_c_in_place_of_N_null": str(wrong_lapse),
            "Hprime_with_Gu3": str(with_center), "wrong_omitted_Gu_Hprime": str(without_center),
            "exact_root_classifications": independent.serialize(classifications),
            "invalid_domain_calls_rejected": rejected}


@cache
def build_report():
    previous = prior_checks()
    residuals = {"literal_lapse_scale_variation_and_Bianchi": potential.checks(),
                 "proper_clocks_and_dynamic_combination": background.checks(),
                 "actual_algebraic_solution_and_original_action": branches.checks()}
    if any(sp.simplify(value) != 0 for group in residuals.values() for value in group.values()):
        raise ValueError("A constant-star background identity failed")
    sources = sorted(ROOT.glob("src/p8_star/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1, "claim": "P8-S6.16.STAR", "date": "2026-09-06",
        "status": "CONSTANT_STAR_SINGLE_NEC_SOURCE_REGULAR_FLAT_NO_BOUNCE; MATCHING_AND_P8_B_OPEN",
        "prior_context_sha256": previous,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
        "formulation": "FORMULATION.md", "written_proof": "notes/proof.md", "source_audit": "notes/sources.md",
        "exact_residuals": {name: dict.fromkeys(group, "0") for name, group in residuals.items()},
        "independent_Fraction_replay": independent.checks(),
        "primary_independent_bridges": independent_bridges(), "checked_controls": controls(),
        "action_and_domain": {
            "potential": "-2 sum_i sqrt|g_i| sum_n beta_in e_n(sqrt(g_i^-1 h)), all constant real beta0..4",
            "kinetics": "EH=-G R_B/2; finite positive leaf G_i, optional constant G_u>=0",
            "curvature": "+---, R_B=-6(DH+2H²), physical Einstein lapse 3G H²=rho",
            "physical_clock": "dT=n_u dt; R_i=a_u/a_i,N_i=n_u/n_i,c_i=R_i/N_i>0; D_tau_i=(R_i/c_i)D_T",
            "geometry": "Finite star, common smooth regular positive spatially flat FLRW/root chart on a connected open interval",
            "source": "Only central homogeneous/isotropic on-shell physical matter; n_h=rho_m+p_m>=0; positive canonical scalar metric sufficient",
            "nontriviality": "G_u>0 or at least one genuine beta1..3 link; endpoint-only leaves are omitted from kinetic combinations",
        },
        "literal_identities": {
            "polynomials": "A=beta0+3beta1R+3beta2R²+beta3R³; B=beta1+3beta2R+3beta3R²+beta4R³; J=beta1+2beta2R+beta3R²",
            "interaction_nulls": "n_leaf=2(R-N)J; n_center=2J(1-c)/R³; n_center+n_leaf/(N R³)=0",
            "unfactored_Bianchi": "6 N_i J_i(R_i)(R_i H_u-H_i)=0; no division by J_i, H_u, null density or auxiliary Hessian",
            "central_null": "G_u H_u'=-n_h/2-sum_i J_i(1-c_i)/R_i³",
            "dynamic_leaf": "R_i'=R_i(1-c_i)H_u; (G_i/R_i²)[H_u'+(1-c_i)H_u²]=J_i(1-c_i)/R_i³",
        },
        "theorem": {
            "statement": "H_u(T0)<=0 implies H_u(T)<=0 for all later T in the same regular interval; no contraction-to-expansion, including degenerate transitions",
            "fixed_dynamic_stratum": "K_D=G_u+sum_dynamic G_i/R_i²>0 gives (H_u/sqrt K_D)'=-n_h/(2 K_D^(3/2))",
            "algebraic_interior": "A genuine J_i has finitely many roots: R_i constant, H_i constant, H_u=c_i H_i/R_i and H_u'=(c_i'/c_i)H_u",
            "arbitrary_zero_sets": "For C1 F_i=H_i-R_i H_u, F_i=F_i'=0 on the complement and hence boundaries of {J_i=0}; no finite-switch premise",
            "compact_comparison": "C=max_genuine sup(|R_i'/R_i|,|c_i'/c_i|)<infinity; H_u'>=0 is not assumed; for H_u>=0, H_u'<=C H_u",
            "global_proof": "Positive part Y=max(H_u,0) is AC and Y'<=CY a.e.; exp(-C(T-T0))Y nonincreasing; exhaust compact intervals",
            "not_global_exact_K": "Algebraic legs cannot generally be included in a single exact all-leg monotonicity identity",
        },
        "actual_controls": {
            "algebraic_solution": "G_u=G_i=1,beta=(0,1,-1/2,0,1/2),a_u=a_i=T^(1/3),n_u=1,n_i=1/(3T),phi=sqrt(2/3)logT,V=0,T>0",
            "algebraic_equations": "A=3/2,B=J=0,H_i=1,H_u=1/(3T); leaf rho,p=3,-3; central interaction zero; scalar rho=p=1/(3T²)",
            "wrong_all_K": "K_all=2 gives 2K_all H_u'-H_u K_all'+n_h=-2/(3T²), not zero",
            "disconnected_exception": "G_u=0,all beta=0,zero matter,Minkowski leaf: arbitrary h and a_u=1+T² with H_u'(0)=2; no physical gravity claim",
            "original_auxiliary_subfamily": "Two leaves beta_i3=p_i,beta_i0=b_i,other beta=0,central -2B det(u),G_u=0; J_i=p_i R_i²",
        },
        "verification_boundary": "Exact symbolic plus separate Fraction/polynomial/first-jet replay and reserved covariant audit support the written topology/positive-part proof; not proof-assistant formalized",
        "primary_sources": ["https://arxiv.org/pdf/1206.3814", "https://arxiv.org/pdf/1211.0214"],
        "not_established": [
            "Any theorem for a disconnected nondynamical G_u=0 center with no genuine link",
            "A globally exact all-leg monotonic quantity across persistent algebraic strata",
            "Scalar/vector/tensor health, a stationary massive spectrum or a Wilsonian cutoff",
            "A statement on nonflat, multimetric-matter, leaf-matter, graph-cycle, derivative-mixing, variable-coefficient or non-NEC models",
            "A continuation through singular, zero or sign-changing lapse/root domains",
            "A UV-positivity/Regge/loop verdict, universal DHOST-row exclusion or original C/D matching",
            "Closure of S6, P8(b), or P8",
        ],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Constant-star no-bounce certificate differs from exact read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.16.STAR: constant-star regular-flat no-bounce replay passed; matching and P8(b) OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
