"""Read-only local-field-domain proportional-vacuum prerequisite replay."""

import argparse
import hashlib
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

import sympy as sp
from p8_variable_beta import verify as prior

from . import bounds, independent, model, obstruction

ROOT = Path(__file__).resolve().parents[2]
P8 = next(path for path in ROOT.parents if path.name == "P8")
REPORT = ROOT/"certificates"/"local-vacuum-obstruction.json"
PRIOR_SHA = "335cd52028baf56b30c75377aa7e6edd7ec86db2799f49b13f467242674fc4d0"
CONTRACT_SHA = "d00df35d5821b05baa61e897725a22ed1ac21d0536b74bdaece3ab0093215901"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serialize(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, sp.Basic):
        if value.has(sp.Float, sp.oo, -sp.oo, sp.zoo, sp.nan):
            raise TypeError("Only exact finite report expressions are permitted")
        return str(value)
    if isinstance(value, dict):
        return {str(key): serialize(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serialize(item) for item in value]
    if value is None or isinstance(value, (str, int, bool)):
        return value
    raise TypeError("Only exact JSON/report values are permitted")


@cache
def prior_checks():
    if sha(prior.REPORT) != PRIOR_SHA or sha(P8/"s6"/"FORMULATION.md") != CONTRACT_SHA:
        raise ValueError("The frozen actual variable-beta action or adopted S6 contract changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {"S6_20_actual_action_and_analytic_local_clock_with_recursive_pins": PRIOR_SHA,
            "adopted_S6_V_and_B_contract": CONTRACT_SHA}


@cache
def independent_bridges():
    records = []
    for fixture in independent.fixtures():
        actual = obstruction.evaluate(*fixture)
        expected = independent.profile(*fixture)
        if {key: Fraction(value) for key, value in actual.items()} != expected:
            raise ValueError("Independent Fraction quotient/clock-derivative replay failed")
        records.append({"u_field_label": str(fixture[0]), "action_c": str(fixture[1]),
                        "independent_ratio_r": str(fixture[2]), "exact_scalar_comparisons": len(expected),
                        "all_equal": True})
    coframes = []
    d = model.derive()
    for fixture in independent.coframe_fixtures():
        ag, ng, r, betas, slopes = fixture
        expected = independent.coframe_fixture(*fixture)
        derivatives = {sp.diff(beta, d["phi"]): sp.Rational(slope) for beta, slope in zip(d["betas"], slopes)}
        at = {d["a_g"]: sp.Rational(ag), d["N_g"]: sp.Rational(ng), d["r"]: sp.Rational(r)}
        at.update({beta: sp.Rational(value) for beta, value in zip(d["betas"], betas)})
        actual = {key: Fraction(d[key].subs(derivatives).subs(at)) for key in expected}
        if actual != expected:
            raise ValueError("Independent literal coframe polynomial variation failed")
        coframes.append({"a_g": str(ag), "N_g": str(ng), "r": str(r),
                         "betas": list(map(str, betas)), "clock_slopes": list(map(str, slopes)),
                         "exact_scalar_comparisons": len(expected), "all_equal": True})
    d, expected = bounds.derive(), independent.coefficients()
    primary = {"polynomial_power_c": tuple(Fraction(sp.Poly(d["Q"], d["c"]).nth(j)) for j in range(3)),
               "upper_at_c2": Fraction(d["Q"].subs(d["c"], 2)),
               "margin_Bernstein": tuple(map(Fraction, d["endpoint_margin_Bernstein"])),
               "c1_U_upper": Fraction(d["constants"]["g_U_upper_c1"]),
               "positive_branch_U_upper": Fraction(d["constants"]["g_U_upper_positive_branch"])}
    if primary != expected:
        raise ValueError("Independently rebuilt continuous polynomial bound failed")
    return {"profile_and_clock_derivative_fixtures": records,
            "literal_coframe_and_clock_variation_fixtures": coframes,
            "independent_continuous_polynomial_and_bound_comparisons": 9,
            "total_exact_scalar_comparisons": sum(row["exact_scalar_comparisons"] for row in records+coframes)+9,
            "independence_boundary": "No SymPy or primary results in the Fraction engine; polynomial margin rebuilt from primitive rational bounds. Finite fixtures check identities, not the continuous theorem by sampling"}


def controls():
    invalid = [lambda args=args: obstruction.evaluate(*args) for args in
               ((0.0, 4), (True, 4), (0, 4.0), (sp.Rational(101, 1000), 4),
                (-sp.Rational(101, 1000), 4), (0, 2), (0, sp.Rational(3, 2)),
                (0, 5), (0, 4, 0), (0, 4, -1), (0, 4, 2.0), (0, 4, True), (sp.oo, 4))]
    invalid += [lambda: independent.profile(0.0, 4), lambda: independent.profile(0, 4, True),
                lambda: serialize(0.1), lambda: serialize(sp.oo)]
    rejected = 0
    for call in invalid:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(invalid):
        raise ValueError("An inexact or out-of-domain input bypassed validation")
    return {"actual_action_omissions_and_changed_action": serialize(obstruction.controls()),
            "generic_equation_system_controls_not_this_action": serialize(model.controls()),
            "invalid_calls_rejected": rejected}


@cache
def build_report():
    previous = prior_checks()
    residuals = {"literal_all_five_beta_metric_and_clock_equations": model.checks(),
                 "actual_local_profile_and_residual_identity": obstruction.checks(),
                 "continuous_rational_bound_identities": bounds.checks()}
    if any(value != 0 for group in residuals.values() for value in group.values()):
        raise ValueError("A local-vacuum exact identity failed")
    audit = ROOT/"tests"/"test_local_vacuum_independent_audit.py"
    if not audit.is_file():
        raise ValueError("The separately authored root audit must exist before freezing")
    sources = sorted(ROOT.glob("src/p8_variable_vacuum/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    d = bounds.derive()
    return {
        "schema": 1, "claim": "P8-S6.28.LOCAL_VACUUM", "date": "2026-09-07",
        "status": "NO_REGULAR_PROPORTIONAL_CONSTANT_CLOCK_MINKOWSKI_VACUUM_IN_SPECIFIED_LOCAL_FIELD_DOMAIN; EXTENSIONS_AND_ORIGINAL_P8_OPEN",
        "prior_context_sha256": previous,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
        "formulation": "FORMULATION.md", "written_proof": "notes/proof.md", "source_audit": "notes/sources.md",
        "exact_symbolic_residuals": {name: dict.fromkeys(group, "0") for name, group in residuals.items()},
        "exact_symbolic_residual_count": sum(map(len, residuals.values())),
        "actual_action_and_domain": {
            "physical_metric": "Same g with positive canonical phi and free chi; both constant in the proposed vacuum, no extra source",
            "convention": "+---; R_B=-6(Hdot+2H²); both Einstein actions -M²R_B/2; interaction -2 sum beta_n e_n",
            "coefficients": "beta_n=M² b_n(u(phi/M),c)/tau², with the actual b2=b3=0",
            "field_interval": "phi=M integral_0^u sqrt(kbar(v,c))dv, |u|<=1/10; strict positive clock gives the analytic inverse",
            "action_parameter": "c=1 or fixed2<c<=4; not the independent proportional vacuum ratio r>0",
            "stationary_geometry": "g=eta,f=r²eta; u is a coefficient label, not vacuum time or a required bounce-to-vacuum history",
        },
        "full_stationary_equations": {
            "g": "U_beta=beta0+3r beta1+3r² beta2+r³ beta3=0; rho_g=2U_beta,pressure_g=-rho_g",
            "f": "V_beta=beta1+3r beta2+3r² beta3+r³ beta4=0; rho_f=2V_beta/r³,pressure_f=-rho_f",
            "clock": "2 sum_n binomial(4,n) r^n beta_n,phi=0, holding r fixed",
            "actual_profile": "b0+3r b1=0; b1+r³b4=0; [2M/(tau²sqrt(kbar))]*(b0,u+4r b1,u+r⁴b4,u)=0",
            "free_chi": "Constant chi automatically solves its unsourced equation and has zero kinetic stress",
        },
        "continuous_proof": {
            "theta": "3h²(c-y)/(2c h'); b4=-(b1/y³)(1-theta)",
            "f_root": "Unique positive root r³=y³/(1-theta); c1 gives r<=y; positive branch gives r<9/4",
            "residual_identity": "U_vac+3y³V_vac/(r²+ry+y²)=3h²/2-nbar/4+9h²y³r³/[2c²(r²+ry+y²)]",
            "constants": serialize(d["constants"]),
            "strict_rational_margins": serialize(d["strict_rational_margins"]),
            "lapse_polynomial": serialize(d["Q"]),
            "lapse_margin_factor": serialize(d["endpoint_margin_factor"]),
            "lapse_margin_Bernstein": serialize(d["endpoint_margin_Bernstein"]),
            "physical_density_gap": "After imposing f-flat exactly: rho_g< -246M²/(25tau²) for c1, and rho_g< -459M²/(400tau²) for2<c<=4",
            "gap_scope": "Not a generic approximate multi-equation error bound; exact f-flat constraint and same local coefficients are essential",
        },
        "independent_replay": independent_bridges(),
        "checked_controls": controls(),
        "adopted_V_B_interpretation": {
            "V": "Algebraic vacuum prerequisite absent within this field domain; no residue, spectrum, amplitude or positivity verdict inferred",
            "B": "A separately named off-interval extension still needs the common-parent domain, state, scales, operators and errors required by the adopted contract",
        },
        "not_established": [
            "No exclusion outside the specified local clock interval; no global real analyticity or unique continuation assumption",
            "Smooth off-interval extensions remain permitted as distinct candidates, not already matched or UV viable",
            "No claim for changed coefficients, additional potentials/sources/heavy-field backgrounds, nonproportional states or other parent actions",
            "No global time trajectory or scattering history connecting bounce and vacuum is required by this argument",
            "No spectrum/positivity/cutoff/loop/UV verdict or original C/D operator match; original P8 remains open",
        ],
        "verification_boundary": "Literal source algebra, independent Fraction/coefficient replay and separately authored audit with continuous written inequalities; not a grid search or proof-assistant formalization",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Local-vacuum certificate differs from exact read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.28.LOCAL_VACUUM: local proportional-vacuum obstruction replay passed; extensions and original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
