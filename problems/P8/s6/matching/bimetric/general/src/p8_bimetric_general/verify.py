"""Read-only arbitrary-beta regular-flat HR no-bounce certificate replay."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_bimetric import verify as prior

from . import background, obstruction, polynomial

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"regular-flat-hr.json"
PRIOR_SHA = "013bc7948f07384ccf6637894c87d73195e5bb8773909f8a8e36c8d41dfda284"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@cache
def prior_checks():
    if sha(prior.REPORT) != PRIOR_SHA:
        raise ValueError("Pinned S6.3.beta1 certificate changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return PRIOR_SHA


def polynomial_bridge():
    bg = background
    variables = {"a": bg.a, "b": bg.b, "Ng": bg.Ng, "Nf": bg.Nf,
                 "ad": sp.diff(bg.a, bg.t), "bd": sp.diff(bg.b, bg.t),
                 "Ngd": sp.diff(bg.Ng, bg.t), "Nfd": sp.diff(bg.Nf, bg.t),
                 "add": sp.diff(bg.a, bg.t, 2), "bdd": sp.diff(bg.b, bg.t, 2),
                 "G": bg.MG2, "F": bg.MF2, "Kg": bg.KG, "Kf": bg.KF,
                 "Vg": bg.VG, "Vf": bg.VF, "m4": bg.M4,
                 **{"beta"+str(index): value for index, value in enumerate(bg.BETAS)}}

    def sympify(value):
        return sum(sp.Rational(coefficient.numerator, coefficient.denominator)
                   *sp.prod(variables[name]**exponent for name, exponent in zip(polynomial.NAMES, power, strict=True) if exponent)
                   for power, coefficient in value.terms.items())

    actual = polynomial.derive()
    return {**{key: sp.factor(sympify(value)-bg.equations()[key]) for key, value in actual["actual"].items()},
            "full_Lagrangian": sp.factor(sympify(actual["L"])-bg.equations()["L"]),
            "full_potential": sp.factor(sympify(actual["potential"])-bg.equations()["interaction"])}


def residuals():
    groups = {"full_two_lapse_and_separate_matter_variation": background.variation_checks(),
              "both_Noether_and_undivided_Bianchi_relations": background.bianchi_checks(),
              "nonzero_P_and_root_cases": obstruction.branch_checks(),
              "physical_CD_window_mismatch": obstruction.window_checks(),
              "explicit_omission_and_outside_domain_points": obstruction.control_checks(),
              "independent_Fraction_to_symbolic_bridge": polynomial_bridge()}
    for name, values in groups.items():
        if any(sp.simplify(value) != 0 for value in values.values()):
            raise ValueError("An arbitrary-beta HR identity failed: "+name)
    return groups


def control_checks():
    values = obstruction.controls()
    zeros = {"zero_lapse_has_zero_metric_determinant", "zero_polynomial_is_a_covered_case_not_a_divisor"}
    for key, value in values.items():
        if (sp.simplify(value) == 0) != (key in zeros):
            raise ValueError("An arbitrary-beta HR control failed: "+key)
    return values


def positive_checks():
    data = obstruction.stationary_equations()
    values = {key: data[key] for key in ("dynamic_positive_inertia", "root_positive_inertia", "second_matter_positive_weight")}
    values.update({"CD_stationary_derivative": 4/obstruction.TAU**2,
                   "CD_window_endpoint_margin": sp.Rational(8, 5)/obstruction.TAU,
                   "CD_window_derivative_margin": sp.Rational(48, 25)/obstruction.TAU**2})
    if any(value.is_positive is not True for value in values.values()):
        raise ValueError("A stipulated-domain positive expression failed")
    return values


@cache
def build_report():
    previous = prior_checks()
    identities, controls, positives = residuals(), control_checks(), positive_checks()
    exact_polynomial = polynomial.exact_checks()
    sources = sorted(ROOT.glob("src/p8_bimetric_general/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {"schema": 1, "claim": "P8-S6.4.HR", "date": "2026-09-06",
            "status": "ARBITRARY_BETA_REGULAR_COMMON_FLAT_HR_NONDEGENERATE_BOUNCE_OBSTRUCTED_WITH_SEPARATE_NEC_MATTER; GENERAL_PARENT_AND_UV_MATCHING_OPEN",
            "prior_S6_3_beta1_sha256": previous,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
            "formulation": "FORMULATION.md", "written_proof": "notes/no-bounce.md",
            "primary_source_audit": "notes/sources.md",
            "primary_sources": {
                "general_FLRW_equations_and_Bianchi_polynomial": "https://arxiv.org/pdf/1111.1655",
                "regular_bimetric_extremality_and_hidden_matter_scope": "https://arxiv.org/pdf/1211.0214",
                "interaction_NEC_anticorrelation_not_physical_matter_NEC": "https://arxiv.org/pdf/1206.3814"},
            "exact_residuals": {name: dict.fromkeys(values, "0") for name, values in identities.items()},
            "independent_Fraction_coefficient_calculation": exact_polynomial,
            "positive_expressions": {key: str(value) for key, value in positives.items()},
            "positive_and_omission_controls": {key: str(value) for key, value in controls.items()},
            "full_algebraic_control_points": {name: {key: str(value) for key, value in values.items()}
                                             for name, values in obstruction.control_points().items()},
            "stationary_equations": {key: str(value) for key, value in obstruction.stationary_equations().items()},
            "theorem_domain": {
                "gravity": "positive G=M_g^2,F=M_f^2; constant arbitrary real beta0,...,beta4; positive m4 (zero interaction included by all beta_i=0)",
                "physical_frame": "g is the fixed physical metric; no added old DHOST action or derivative metric-frame change",
                "geometry": "both metrics have a common spatially flat FLRW foliation and positive-root square-root eigenvalues (c,y,y,y)",
                "regularity": "a,b,Ng,Nf smooth and finite, positive and nonzero in a neighborhood of the stationary slice",
                "primary_matter": "all minimally coupled classical NEC matter on g, no matter on f",
                "separate_sector_extension": "two independently minimally coupled, separately conserved matter sectors, NEC in their own metrics; no shared doubly coupled field",
                "target": "H_g=0 with D_g H_g>0, including CD; no vacuum/mass-gap/perturbation-health assumption for arbitrary beta"},
            "full_case_split": {
                "P_nonzero": "P(y0)!=0 gives a neighborhood with P!=0 by continuity; undivided Bianchi then yields Hf=Hg/y there; at stationarity DfHf=Hdot_g/(c*y)",
                "P_zero": "at any P(y0)=0, g interaction null stress vanishes and -2G Hdot_g=rho_g+p_g>=0 directly; never divide by P or infer a dynamic branch relation",
                "coverage": "simple/double/isolated roots, branch-switching points, algebraic intervals and the identically zero polynomial are included pointwise"},
            "nonzero_P_weighted_equation": "-2*(G+F*y^2)*Hdot_g=(rho_g+p_g)+c*y^3*(rho_f+p_f)",
            "root_equation": "-2*G*Hdot_g=rho_g+p_g",
            "robust_CD_mismatch": "On [-tau/2,tau/2] in physical g cosmic time, simultaneous endpoint H errors <8/(5tau) and uniform Hdot error <48/(25tau^2) are impossible",
            "scope_of_controls": "forbidden-NEC jets are algebraic controls, not healthy solutions; formal c<0 continuation only tests positivity failure, not the original positive-root action",
            "verification_boundary": [
                "Pinned beta1 and both of its input lineages are replayed without prior-file edits",
                "Arbitrary-coupling and jet identities are coefficientwise exact, not rational sampling or interval interpolation",
                "Independent full covariant/source/branch audit is hashed and tested separately",
                "Pointwise root-case coverage, continuity and intermediate-value implications have written proofs, not formalized global PDE proofs",
                "Literature sources supply compatible identities and scope cautions; no global novelty claim is made"],
            "not_established": [
                "Any healthy arbitrary-beta Minkowski or rolling spectrum; arbitrary-beta vacuum matching is not inherited from beta1",
                "An exclusion of non-bidiagonal/non-common-FLRW metrics, non-flat spatial curvature, singular or other square-root branches",
                "An exclusion of quantum/NEC-violating/nonminimal/derivative/shared doubly coupled matter or other metric-frame proposals",
                "An exclusion of every degenerate bounce, all bimetric modifications or general UV completions",
                "A controlled alternative parent, radiative closure, finite-gravity positivity, or a physical interaction cutoff",
                "Completion of S6, P8(b), P8(a), or P8"]}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Arbitrary-beta HR certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.4.HR: arbitrary-beta regular-flat nondegenerate-bounce obstruction replay passed, including roots; general matching OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
