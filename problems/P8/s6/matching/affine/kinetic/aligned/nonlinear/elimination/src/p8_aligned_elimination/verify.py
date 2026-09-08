"""Read-only ordered vector elimination and full nonlinear source remainder."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_affine_nonlinear import verify as parent

from . import bounds, operators, source

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"ordered-elimination-source-remainder.json"
PARENT_SHA = "11c833b6cad1c707da69de3a6852eaeba6cbb20e23175d712e7374ac572f0b63"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_aligned_elimination/*.py"))
            +sorted(ROOT.glob("tests/*.py"))+sorted(ROOT.glob("*.md"))
            +sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen local nonlinear-constraint certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_45_fully_rebuilt": PARENT_SHA,
            "actual_full_action_ODE_source_prepared_light_and_all64_ancestry_rechecked": True,
            "no_frozen_action_or_ancestor_modified": True}


def residuals():
    out = {**source.checks(), **operators.checks()}
    out["prepared_full_source_jet_through_three"] = source.preparation()["prepared_full_source_jet_through_three"]
    return out


def controls():
    bad = (True, False, 1.0, sp.Float(1), "1", sp.I, sp.oo, sp.nan, sp.Symbol("unproved"))
    calls = [lambda value=value: bounds.values(value, 1) for value in bad]
    calls += [lambda value=value: bounds.values(sp.Rational(1, 100), value) for value in bad]
    calls += [lambda value=value: bounds.values(value, 1) for value in (0, -1, sp.Rational(1, 20))]
    calls += [lambda value=value: bounds.values(sp.Rational(1, 100), value) for value in (0, -1)]
    calls += [lambda value=value: bounds.physical_values(sp.Rational(1, 100), 1, 1, value) for value in bad]
    for name in ("mass_squared", "time_scale"):
        calls += [lambda value=value, name=name: bounds.physical_values(sp.Rational(1, 100), 1, **{name: value})
                  for value in (0, -1)]
    calls += [lambda value=value: bounds.physical_values(sp.Rational(1, 100), value) for value in (0, -1)]
    calls += [lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan)]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An inexact or outside-source-domain input was accepted")
    return {"rejected_inputs": rejected,
            "noncommuting_mass_and_curl_operators_kept": True,
            "one_copy_retarded_action_countercontrol": True,
            "full_source_is_not_a_finite_Fourier_polynomial": True,
            "Q_remainder_spatial_chain_rule_not_omitted": True,
            "source_bound_not_an_inverse_or_nonlinear_solution_bound": True,
            "quantum_determinant_identified_not_evaluated_or_bounded": True}


@cache
def build_report():
    pins = prior_checks()
    identities = residuals()
    exact = affine.certify_residuals(identities)
    proofs = {**bounds.proof_checks(), **operators.proof_checks()}
    if not all(value is True for value in proofs.values()):
        raise ValueError("An ordered, continuous source-bound or causal-interface check failed")
    return {"schema": 1, "claim": "P8-S6.46.ALIGNED_ELIMINATION", "date": "2026-09-07",
            "status": "EXACT_ORDERED_VECTOR_ELIMINATION_AND_NONLINEAR_SOURCE_REMAINDER_CERTIFIED; ORIGINAL_P8_OPEN",
            "prior_sha256": pins,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md", "written_proofs": ["notes/source.md", "notes/operators.md"],
            "exact_residuals": exact, "named_exact_check_count": len(identities),
            "checked_scalar_entries": sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1
                                          for value in identities.values()),
            "proof_checks": proofs,
            "literal_action": "Unchanged S6.42; at fixed light fields the full retained vector action is 1/2<W-S,M(W-S)>+zeta/2<W,DW>, with actual original ODE source and physical metric",
            "ordered_elimination": "M-M*K^-1*M=zeta*D-zeta^2*D*K^-1*D, K=M+zeta*D; exact common-domain identity without commuting M,D or asserting an inverse norm",
            "state_and_boundary": "A stationary induced action uses a compatible symmetric variational domain. Retarded response keeps its homogeneous data and constraint compatibility; one-copy variation of G_ret symmetrizes it and introduces an advanced contribution.",
            "full_source_domain": "All real u, r=N-1 and k=K_hat-3H; |r|,|k| and physical spatial gradient norms <=epsilon<=1/25, normalized with tau. This ensures the original X tube, not a quantitative S6.45 nonlinear neighborhood or full on-shell solution.",
            "actual_ODE_remainders": "|Q-f1*(x+1)^2/2|<=(140/93)|x+1|^3/h^2; |Q_x-f1*(x+1)|<=(140/31)|x+1|^2/h^2, f1=-9H/(8h^2)",
            "continuous_constants": serialize(bounds.constants()),
            "source_and_local_density_bounds": "Normal remainder 104epsilon^3, coordinate remainder 75epsilon^3, coordinate spatial-gradient remainder 230epsilon^3, electric remainder 300epsilon^3; same-physical-metric local curl-density substitution error <=14325*zeta*epsilon^5",
            "preparation": "r=r'=r''=0 gives S=S'=S''=0 and S'''=-2r'''k_initial/h; the nontrivial S6.44 profiles retain source preparation as specified fields, not as a claimed full nonlinear solution",
            "Fourier_boundary": "The full rational lapse map produces further harmonics and variable vector coefficients; the old diagonal fixed-output background response estimate is not extended here",
            "normalized_example": serialize(bounds.values(sp.Rational(1, 1000), sp.Rational(1, 10**9))),
            "physical_example": serialize(bounds.physical_values(sp.Rational(1, 1000), sp.Rational(1, 10000), 3, 2)),
            "quantum_boundary": "The regulated vector Gaussian leaves det(K)^(-1/2); source alignment does not remove its light-dependent determinant. No measure/complement/light loop, subtraction, finite counterterm or quantum bound is supplied.",
            "controls": controls(),
            "verdict": "The exact classical induced operator and genuine full-source/spatial-curl remainders are established. The nonlocal inverse, causal nonlinear feedback and quantum corrections remain distinct unresolved estimates.",
            "not_established": ["A full nonlinear solution, variable-coefficient inverse or nonlocal induced-action/variation remainder bound",
                                "Nonlinear health, a quantified S6.45 neighborhood, quantum corrections, interacting cutoff or stationary gap",
                                "V/G/B UV admissibility or original P8 closure"],
            "verification_boundary": "Exact operator ordering, independent dense stationary solves, causal countercontrol and written continuous original-ODE/source inequalities; not proof-assistant formalization"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The ordered elimination/source remainder report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.46.ALIGNED_ELIMINATION replay passed; exact identity and source bounds only, original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
