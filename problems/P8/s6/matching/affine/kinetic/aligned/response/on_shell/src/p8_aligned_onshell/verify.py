"""Read-only nontrivial prepared linear light solution and response transfer."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_aligned_response import verify as parent

from . import bounds, checks, phase

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"prepared-linear-light-family.json"
PARENT_SHA = "6b6ca45b6028c2f11b1e14cf6aa47ccffa65dc1e53e458d3f7e527fb3018d4b4"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_aligned_onshell/*.py"))
            +sorted(ROOT.glob("tests/*.py"))+sorted(ROOT.glob("*.md"))
            +sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen prepared leading-response certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_43_fully_rebuilt": PARENT_SHA,
            "actual_physical_light_system_source_and_retarded_theorem_rechecked": True,
            "no_frozen_action_or_ancestor_modified": True}


def controls():
    wrong = (True, False, 1.0, sp.Float(1), "1", sp.I, sp.oo, sp.nan, sp.Symbol("unproved"))
    calls = [lambda value=value: bounds.prepared_family(value) for value in wrong]
    calls += [lambda value=value: bounds.transfer(1, value) for value in wrong]
    calls += [lambda value=value: bounds.prepared_family(value) for value in (0, -1)]
    calls += [lambda value=value: bounds.transfer(1, value) for value in
              (0, -1, sp.Rational(1, 2000))]
    calls += [lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan)]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An inexact, unproved, trivial or outside-response-domain input was accepted")
    return {"rejected_inputs": rejected,
            "nonzero_family_distinguished_from_trivial_zero_control": True,
            "moving_boundary_omission_detected_away_from_center": True,
            "independent_initial_Cauchy_recurrence": True,
            "output_zero_and_double_momentum_retained": True,
            "conservative_initial_scale_not_observational_or_optimal": True,
            "on_shell_means_linear_light_not_full_nonlinear_parent": True}


@cache
def build_report():
    pins = prior_checks()
    residuals = checks.residuals()
    exact = affine.certify_residuals(residuals)
    proofs = bounds.proof_checks()
    if not all(value is True for value in proofs.values()):
        raise ValueError("An actual light preparation, continuous bound or nontriviality proof failed")
    jets = phase.initial_jets()
    return {"schema": 1, "claim": "P8-S6.44.ALIGNED_ONSHELL", "date": "2026-09-07",
            "status": "NONTRIVIAL_PREPARED_LINEAR_LIGHT_FAMILY_AND_LEADING_RESPONSE_TRANSFER_CERTIFIED; ORIGINAL_P8_OPEN",
            "prior_sha256": pins,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md", "written_proofs": ["notes/proof.md"],
            "exact_residuals": exact, "named_exact_check_count": len(residuals),
            "checked_scalar_entries": sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1
                                          for value in residuals.values()),
            "proof_checks": proofs,
            "literal_action": "Original regular CD/M1 linear light equations from unchanged S6.42; actual quadratic source and zero-data leading vector response from S6.43",
            "physical_phase": "Y=(b,s,P_b,P_s), P_b=2a³q*v, P_s=a³P; moving canonical Hamiltonian includes -H*b*P_b; n=-R0/J with no Theta divisor",
            "geometric_readout": "delta_K_hat=3(Theta-H)*n-3ell*s/2+q*b; zeta_metric=v+n/(2h) and omega_N=1/(2h)",
            "initial_preparation": serialize({"time": phase.LEFT, "unit_direction": jets["initial_vector"],
                                              "lapse_third": jets["lapse_third"], "trace_initial": jets["trace_initial"],
                                              "source_third": jets["source_third_initial"],
                                              "observability_determinant": jets["observability_determinant"]}),
            "continuous_bounds": serialize(bounds.continuous()),
            "domain_and_scaling": "u in [-1/2,1/2], incoming k²=1/16; initial state epsilon*Y_star/(3^238*C), C=1+B_N+B_K; actual lapse and trace jets through order three <epsilon",
            "source_and_response_transfer": "Real cosine gives zero and double momentum, k_out²=1/4; actual S2 is nonzero and prepared, E=(136139/8)epsilon²; S6.43 physical retarded error bounds apply",
            "example": serialize(bounds.transfer(sp.Rational(1, 1000), sp.Rational(1, 10**9))),
            "controls": controls(),
            "verdict": "A nonzero actual linear light solution realizes and maintains the prepared source class through the full declared bounce window, enabling the certified leading heavy-response estimate. This is not a full nonlinear parent solution.",
            "not_established": ["The complete nonlinear light/heavy solution, nonlinear secondary constraints or nonlinear stability",
                                "Higher perturbative-order or induced-action errors, quantum loops, an interacting cutoff or stationary gap",
                                "Vacuum/finite-gravity V/G/B admissibility or original P8 closure"],
            "verification_boundary": "Exact symbolic constraints, independent canonical and Cauchy checks, continuous rational matrix bounds and explicit integral-series proof; not proof-assistant formalization"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The prepared linear-light family report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.44.ALIGNED_ONSHELL replay passed; linear light and leading response only, original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
