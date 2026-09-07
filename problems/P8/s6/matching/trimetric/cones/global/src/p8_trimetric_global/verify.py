"""Read-only exact replay of the actual auxiliary-background no-bounce gate."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_trimetric_cones import verify as prior

from . import background, guards, independent, monotonic, obstruction, vacuum

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"auxiliary-no-bounce.json"
P8 = next(path for path in ROOT.parents if path.name == "P8")
PRIOR_SHA = "ab280acf8d20fb2ae18d75c265deecd59832e993ba0338e9e5c8f7691c0a7e6e"
P8_FORMULATION_SHA = "73dd18f7a4b56a0205f9fbc6a4b09b213a3c280b3f0274414a255635b04e2b09"
CD_SHA = "caf8c8e688a7565b9d00f921c099a28da00f97522ed26ad182a8227eb80cd4dd"
S6_FORMULATION_SHA = "d00df35d5821b05baa61e897725a22ed1ac21d0536b74bdaece3ab0093215901"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@cache
def prior_checks():
    pins = {"S6_14_actual_auxiliary_action_and_cones": (prior.REPORT, PRIOR_SHA),
            "original_B_curvature_and_physical_metric": (P8/"FORMULATION.md", P8_FORMULATION_SHA),
            "conditional_original_CD_M1_endpoint_target": (P8/"certificates"/"witness-CD_matter.json", CD_SHA),
            "adopted_S6_controlled_matching_contract": (P8/"s6"/"FORMULATION.md", S6_FORMULATION_SHA)}
    if any(sha(path) != expected for path, expected in pins.values()):
        raise ValueError("A pinned action, physical dictionary or matching contract changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {key: expected for key, (_, expected) in pins.items()}


def plain(value):
    if isinstance(value, dict):
        return {key: plain(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [plain(item) for item in value]
    return str(value) if isinstance(value, sp.Basic) else value


def controls():
    v, g, m, o = vacuum.controls(), guards.controls(), monotonic.controls(), obstruction.controls()
    if tuple(v.values()) != (1, -2, -1, 2, 3, 4, sp.Rational(10, 3)):
        raise ValueError("The non-unit actual-vacuum mass dictionary failed")
    t = sp.Symbol("T", real=True)
    if tuple(g.values()) != (-2, 0, 0, sp.Rational(3, 2), sp.Rational(37, 12), 0, 2, sp.Rational(1, 2), -6*t**2):
        raise ValueError("An actual-versus-algebraic or disconnected/singular guard failed")
    d = monotonic.derive()
    if sp.simplify(m["positive_bounce_H0_residual"]-2*d["K"]-d["n_h"]) != 0:
        raise ValueError("The exact regular positive-bounce residual failed")
    if tuple(m.values())[1:] != (3, sp.Rational(5, 2), -2, -6*t**2, sp.Rational(8, 5), 2):
        raise ValueError("A residual-budget or actual CD endpoint threshold failed")
    if [entry["strict_chord_margin"] for entry in o.values()] != [1, 0, -sp.Rational(3, 5), sp.Rational(7, 10), 0]:
        raise ValueError("A finite negative-link three-slice control failed")
    if monotonic.cd_endpoint_test(sp.Rational(1, 2), sp.Rational(8, 5))["status"] != "INCONCLUSIVE":
        raise ValueError("An endpoint equality was treated as a strict obstruction")
    for bad in (True, 0.1, sp.Float("0.1"), sp.oo, sp.I):
        try:
            monotonic.cd_endpoint_test(bad, 0)
        except TypeError:
            pass
        else:
            raise ValueError("An inexact or invalid certificate-bound input was accepted")
    return plain({"nonunit_physical_vacuum_mass": v, "actual_solution_and_domain_guards": g,
                  "monotonicity_and_actual_endpoint_bounds": m, "complementary_three_slice_tests": o})


def dynamics_bridge():
    """Primary expressions versus independent literal Dual differentiation."""
    result = []
    for fixture in independent.phase_fixtures():
        d = monotonic.derive(fixture["links"])
        inputs = fixture["inputs"]
        point = {d[key]: sp.Rational(inputs[key]) for key in ("H_u", "H_u_prime", "n_h")}
        for link, values in zip(d["links"], inputs["links"], strict=True):
            point.update({link[key]: sp.Rational(value) for key, value in values.items()})
        outputs = {"K": "K", "K_prime": "K_prime", "combined_residual": "combined_residual",
                   "normalized_Hubble_prime": "Z_prime"}
        for primary_key, independent_key in outputs.items():
            if sp.simplify(d[primary_key].subs(point)-sp.Rational(fixture[independent_key])) != 0:
                raise ValueError("Primary versus independent source/dynamics/normalization fixture failed")
        result.append({"inputs": inputs, "primary_minus_independent": dict.fromkeys(outputs, "0")})
    return result


def rolling_bridge():
    """Compare actual positive-NEC local solution jets, including q<0."""
    d = monotonic.derive()
    result = []
    for fixture in independent.actual_rolling_fixtures():
        values = {key: sp.Rational(value) for key, value in fixture.items()}
        h = sp.sqrt(values["H_u_squared"])
        point = {d["H_u"]: h, d["H_u_prime"]: values["H_u_prime"], d["n_h"]: 2*values["rho"]}
        for link in d["links"]:
            point.update({link["G"]: 1, link["R"]: values["A"], link["c"]: values["A"]/values["N"],
                          link["p"]: values["q"]})
        residuals = {"K": d["K"].subs(point)-values["K"],
                     "K_prime": d["K_prime"].subs(point)-h*values["K_prime_over_H_u"],
                     "combined_residual": d["combined_residual"].subs(point),
                     "normalized_equation_residual": d["normalized_equation_residual"].subs(point)}
        if any(sp.simplify(value) != 0 for value in residuals.values()):
            raise ValueError("The primary equation failed an independently reconstructed actual local solution")
        result.append({"A": str(values["A"]), "q": str(values["q"]),
                       "primary_minus_independent": dict.fromkeys(residuals, "0")})
    return result


@cache
def build_report():
    previous = prior_checks()
    residuals = {"literal_full_covariant_and_lapse_retaining_background": background.checks(),
                 "two_nonzero_link_Bianchi_and_physical_monotonicity": monotonic.checks(2),
                 "one_nonzero_link_Bianchi_and_physical_monotonicity": monotonic.checks(1),
                 "conditional_flat_vacuum_and_full_FP_invariant": vacuum.checks(),
                 "complementary_cone_cap_and_three_slice_identities": obstruction.checks(),
                 "actual_singular_disconnected_and_deSitter_controls": guards.checks()}
    if any(sp.simplify(value) != 0 for group in residuals.values() for value in group.values()):
        raise ValueError("A full-equation, monotonicity, vacuum or scope identity failed")
    for count in (1, 2):
        d = monotonic.derive(count)
        if d["K"].is_positive is not True or any(link["p"] in d["K"].free_symbols for link in d["links"]):
            raise ValueError("Positive normalization must not depend on link signs or a TT inverse")
    rational, omissions = independent.checks(), controls()
    sources = sorted(ROOT.glob("src/p8_trimetric_global/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    if ROOT/"tests"/"test_auxiliary_global_independent_audit.py" not in sources:
        raise ValueError("The separately authored full covariant audit is required")
    return {
        "schema": 1, "claim": "P8-S6.15.AUXILIARY", "date": "2026-09-06",
        "status": "ALL_LINK_SIGNS_ACTUAL_PHYSICAL_HUBBLE_MONOTONICITY; SPECIFIED_AUXILIARY_PARENT_NO_BOUNCE; P8_OPEN",
        "prior_context_sha256": previous,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
        "formulation": "FORMULATION.md", "written_proof": "notes/proof.md", "source_audit": "notes/sources.md",
        "exact_residuals": {name: dict.fromkeys(group, "0") for name, group in residuals.items()},
        "independent_Fraction_replay": rational,
        "independent_source_dynamics_normalization_bridge": dynamics_bridge(),
        "independent_actual_local_solution_bridge": rolling_bridge(),
        "checked_controls": omissions,
        "action_and_domain": {
            "action": "Sum -Gi/2 det(ei) R[gi] -2det(u)[B+tr(u^-1 sum pi ei)] -2sum bi det(ei) +epsilon Sm[h], h=u^T eta u",
            "constants": "Gi>0; pi,bi,B constant real of either sign; I={i:pi!=0} nonempty; one or two active links, no disconnected metric in K",
            "excluded_operators": "No u Einstein/other derivative terms, other pairwise potential, direct e/v matter, or added derivative source operators",
            "B_curvature": "+---; R_B=-6(DH+2H²), Einstein_B00=+3H², EH=-Gi R_B/2; not the opposite A/FK convention",
            "actual_geometry": "Connected smooth regular common spatially flat homogeneous interval, positive lapse/scale entries and the parent coframe/root domain",
            "physical_clock": "T is actual h proper time; Ri=au/ai,ci=Ri ni/nu,d taui/dT=ci/Ri>0",
            "physical_source": "Actual u-only isotropic n_h=epsilon(rho+p)>=0; positive canonical matter is sufficient but strict rolling is unnecessary",
            "conservation": "Full Einstein equations imply the interaction Bianchi factor, not an extra independent interaction NEC or twice-assumed matter conservation",
            "units": "Gi,K mass²; pi,bi,B,n_h mass4; H_u mass; Kprime mass3; H_u/sqrtK dimensionless; no implicit old-M1 Planck dictionary",
        },
        "primary_full_background_theorem": {
            "full_Einstein_equations": "3Gi Hi²=2pi Ri³+2bi; Gi(2DiHi+3Hi²)=2pi Ri³/ci+2bi",
            "undivided_source_Bianchi": "Di rhoi+3Hi(rhoi+pressure_i)=6pi Ri³(Ri Hu-Hi)/ci=0, hence Hi=Ri Hu for every active link, even at Hi=Hu=0",
            "ratio_and_positive_normalizer": "Ri'=Ri(1-ci)Hu; K=sum_active Gi/Ri²>0; K'=-2Hu sum_active Ki(1-ci)",
            "source_residual": "U=sum_active(pi/Ri)(ci-1)-n_h/2",
            "physical_Einstein_residuals": "Ei=Ki[Hu'+(1-ci)Hu²]-(pi/Ri)(1-ci)",
            "exact_residual_combination": "2sum Ei-2U=2K Hu'-Hu K'+n_h",
            "on_shell_identity": "2K Hu'-Hu K'=-n_h; (Hu/sqrtK)'=-n_h/(2K^(3/2))<=0",
            "conclusion": "No negative-to-positive physical Hubble sign change on the connected regular interval, including every degenerate bounce",
            "no_extra_premises": "No cone condition, flat vacuum, tail/completeness assumption, Hubble division, uniform tail bound on K, or auxiliary TT inverse",
            "singular_TT_scope": "Includes sum pi/Ri=0 background points/intervals; no singular TT inverse is used or claimed",
        },
        "conditional_actual_CD_matching": {
            "required_dictionary": "An actual full-parent solution and matching identification of its actual h/proper T with the original CD target at both endpoints; reduced approximate equations alone do not supply this",
            "target": "aCD=(1+(T/tau)²)², HCD=4T/(tau²+T²),tau>0",
            "necessary_maximum_endpoint_error": "max_± |Hu(±Ltau)-HCD(±Ltau)| >=4L/[tau(1+L²)],L>0",
            "half_duration": "8/(5tau)", "unit_duration": "2/tau", "strictness": "Equality is inconclusive; no Hubble-derivative error premise is required",
            "light_only_scope": "Constrains even an asymmetric light-only proposal if it claims this controlled actual full-parent background match; not a universal reduction or source/operator matching theorem",
            "matter_and_vacuum_guards": "No identification of one internal scalar with old DHOST clock plus free chi, and no physical vacuum-to-bounce trajectory premise",
        },
        "conditional_actual_residual_budget": {
            "hypotheses": "K>=kappa>0,Hu'>=a>0,n_h>=nu>=0,|Hu|<=eta,|K'|<=D",
            "lower_bound": "2K Hu'-Hu K'+n_h>=2kappa*a+nu-eta*D",
            "meaning": "A positive bound is a necessary combined actual-equation correction cancellation, not a computed omitted-operator, loop, field-map or cutoff error estimate",
            "numeric_domain": "Exact finite rational values and explicit sign gates; reject bool, ordinary/SymPy Float, nonfinite and nonreal values",
        },
        "complementary_individual_cone_tail_proof": {
            "extra_contract": "Every active individual Einstein cone ci<=1; not a hidden premise of the primary theorem",
            "negative_link": "pi<0 implies bi>0 and Ri<=Rcap=(bi/|pi|)^(1/3), while DiHi<=0; bi<=0 has no regular real flat-FLRW point",
            "three_slice": "log ai is concave in ordered proper taui, ai>=au/Rcap; two finite endpoints above any fixed middle value contradict the chord inequality",
            "global_tails": "If au grows without bound at both physical tails, the finite endpoints exist; no infinite Einstein proper-time premise and no vacuum needed",
            "finite_CD_test": "r=Ri(0)/Rcap in(0,1],L>0; r(1+L²)²(1-epsa)>1+epsa,0<=epsa<1, is a sufficient strict three-slice obstruction",
            "nonnegative_links": "u null identity +NEC+ci<=1 force n_h=0 and ci=1 on active links; constant Hi and constraint fix constant Ri, so Hu is constant",
        },
        "separate_same_action_vacuum_and_TT_inverse": {
            "premise": "Actual positive Lorentz-flat proportional vacuum of this same action, zero scalar gradient and full u/scalar equations, then ei0=u0/Ri0 and bi=-pi Ri0³",
            "constant_potential": "Beffective=B+epsilon V0/2=-3S0,Pi0=pi/Ri0,S0=sumPi0; endpoint calibration or positive cap alone does not prove a vacuum",
            "literal_stationary_density": "-2S0 det(I+(Pg0 H+Pf0 J)/(2S0))+2Pg0 det(I+H/2)+2Pf0 det(I+J/2),S0!=0",
            "full_traceful_FP_quadratic": "Pg0 Pf0/(4S0)*[(tr(H-J))²-tr((H-J)²)]",
            "actual_physical_vacuum_mass": "Ki0=Gi/Ri0²>0,qeff0=2Pg0Pf0/S0,mFP0²=qeff0*(1/Kg0+1/Kf0)",
            "sign_free_weight_loss": "Pi-Pi0=-3Gi Hi²/[2Ri Ri0(Ri²+Ri Ri0+Ri0²)]<=0 from the calibrated lapse residual",
            "mixed_positive_FP": "Genuine mixed links and finite positive mFP0² imply S0<0, hence S<=S0<0 throughout every regular same-action flat-FLRW point",
            "physical_full_tensor_corollary": "For homogeneous canonical scalars without independent TT matter response: KTT_i=Gi/(ci Ri²),FTT_i=Gi ci/Ri²; nonsingular algebraic TT elimination adds no derivatives and full squared cones are ci²",
            "limits": "No scalar/vector health, positive rolling Green inverse, gap/adiabatic control, light-only or universal vacuum inference; primary no-bounce needs none of these additional premises",
        },
        "nonempty_and_outside_hypothesis_controls": {
            "actual_singular_flat": "pg=1,pf=-1,B=0,bg=-1,bf=1,e=v=u=I,zero source solves full equations but S=0 and uTT imposes gamma_g=gamma_f; main background theorem still applies",
            "actual_all_zero_links": "All links/endpoints/B/V=0,constant free scalar,e=v=I,u=(1+t²)I solves all equations; T=t+t³/3 and Hu'(0)=2, but I empty and the physical auxiliary metric is undetermined",
            "actual_deSitter": "Gi=pi=1,B=-6,bi=-5/8,e=v=u with H=1/2,constant free scalar solves all equations and demonstrates allowed zero-NEC expansion",
            "actual_local_rolling_family": "Gi=epsilon=1,pi=q,B=-6q,bi=-q,V=0,e=v=r; N=A/(6A-5),rho=p=12q(A-1)/A,Hr²=2q(A³-1)/3,dA/dt=-6Hr A(A-1)/(6A-5)",
            "rolling_domain": "Positive Hr root,a_r>0,A>5/6,q(A-1)>0; q<0 strictly 5/6<A<1 and q>0 requires A>1; analytic local ODE, no uniform duration or healthy negative-q vacuum-spectrum claim",
            "rolling_full_equations": "a_r'=Hr a_r,psi'=N sqrt(2rho),Hr'=-6qA³(A-1)/(6A-5),Hu=Hr/A; all lapses, accelerations, u equations, scalar current and normalized identity hold",
            "calibration_only_not_vacuum": "pg=-2,pf=1,bg=2,bf=-1,B=4 at identity coframes satisfies Einstein calibration but fails the actual u equation",
            "degenerate_kinematic_only": "K=1,Hu=T³ requires n_h=-6T²; this is an exact NEC diagnostic, not a full parent solution",
        },
        "primary_source_audit": {
            "1804_04671": {"url": "https://arxiv.org/pdf/1804.04671", "use": "Chosen auxiliary action/source equations and actual-vierbein distinction; endpoint extension explicit; no cutoff imported"},
            "1211_0214": {"url": "https://arxiv.org/pdf/1211.0214", "use": "Prior proper-velocity Bianchi lock and extremum relations, eq15/18-19; differing branch restrictions and arbitrary hidden matter, not our u-only NEC theorem"},
            "1206_3814": {"url": "https://arxiv.org/pdf/1206.3814", "use": "Section4 interaction effective-NEC anticorrelation, explicitly distinct from actual physical matter NEC"},
            "novelty": "No general novelty claim; scoped full-source theorem derived here, not imported from prior papers or a mutable successor",
        },
        "verification_boundary": "Exact symbolic plus independent coefficientwise Fraction/Dual checks and genuine local-solution fixtures, separately authored covariant action audit and written global calculus proof; not proof-assistant formalized or finite sampling of all time",
        "not_established": [
            "A no-bounce theorem for extra potential/derivative operators, u kinetics, direct e/v matter, non-NEC/quantum sources, anisotropic/curved or singular geometries",
            "The theorem for an empty nonzero-link index set or an invertible auxiliary construction in that underdetermined case",
            "A vacuum from endpoint calibration alone or a regular TT inverse at S=0 without the separate calibrated positive-FP premise",
            "A scalar/vector stability, positive heavy propagator, finite-band light-only control, cutoff or loop theorem from the vacuum mass or full tensor cones",
            "A controlled actual CD matching dictionary or full matter/operator match merely from reduced equations",
            "A universal exclusion across every parent or a replacement of adopted S6 by full UV completion",
            "An alteration to the completed scoped photon objective or frozen linear classification; closure of adopted S6 or original P8",
        ],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Actual auxiliary-background no-bounce certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.15.AUXILIARY: all-sign actual physical-Hubble no-bounce replay passed; original matching and P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
