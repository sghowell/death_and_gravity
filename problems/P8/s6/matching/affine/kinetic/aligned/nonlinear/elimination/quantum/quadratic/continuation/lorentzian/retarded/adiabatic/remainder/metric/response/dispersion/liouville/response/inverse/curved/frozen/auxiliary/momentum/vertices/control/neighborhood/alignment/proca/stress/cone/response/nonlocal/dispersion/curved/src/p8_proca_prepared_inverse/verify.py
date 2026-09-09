"""Read-only actual prepared ordinary-Proca tree-plus-loop inverse."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as original
from p8_proca_rank_one_inverse import verify as parent

from . import audit, coordinates, matching, tree, volterra

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"prepared-ordinary-Proca-inverse.json"
PARENT_SHA="0f04fd277040ad21525eb1d9102e6af45c4d695fac3f69d0d37710a7ab8a0db6"
sha,serialize=original.sha,original.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_proca_prepared_inverse/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen rank-one ordinary Proca scalar inverse report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_86_fully_rebuilt":PARENT_SHA,
            "actual_new_curved_currents_tree_and_fixed_profile_rederived":True}


@cache
def build_report():
    prior=prior_checks()
    identities,gates=audit.residuals(),audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A new prepared ordinary Proca coupled inverse gate failed")
    exact=original.certify_residuals(identities)
    return {"schema":1,"claim":"P8-S6.87.PREPARED_ORDINARY_PROCA_INVERSE","date":"2026-09-09",
            "status":"ACTUAL_NEW_PREPARED_HOMOGENEOUS_TREE_PLUS_FIXED_PROFILE_ONE_LOOP_RESPONSE_HAS_UNIQUE_CAUSAL_INVERSE_WITH_CONTINUOUS_ORIGINAL_METRIC_RECONSTRUCTION; QUANTUM_STABILITY_MATCHING_AND_ORIGINAL_P8_OPEN",
            "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
            "formulation":"FORMULATION.md",
            "written_proofs":["notes/coordinates.md","notes/tree.md","notes/matching.md","notes/volterra.md","notes/scope.md"],
            "exact_residuals":exact,"named_exact_check_count":len(identities),
            "checked_scalar_entries":sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in identities.values()),
            "proof_checks":gates,
            "actual_new_time_covariance_and_frozen_profile_currents":serialize(coordinates.data()),
            "actual_literal_tree_zero_charge_reduction_and_all_derivative_coefficients":serialize(tree.data()),
            "new_full_dimensional_current_pairs":serialize({sector:matching.high(sector) for sector in ("T","L")}),
            "new_curved_scalar_leading_log_and_fixed_finite_contact_matching":serialize(matching.data()),
            "actual_prepared_block_inverse_and_original_force_reconstruction":serialize(volterra.data()),
            "exact_fixed_mass_coupling_normalizations":serialize({value:volterra.coupling(value) for value in (1,2,400)}),
            "controls":audit.controls(),
            "verdict":"The actual new selected-state response has one nonlocal scale channel after an initially identity time change, while the fixed scalar profile supplies retained local time contacts. The unchanged literal scalar-matter action supplies a uniformly nonzero fourth-order time channel. New dimensional current and finite-local matching gives I4(T+gamma Q)=diag(-6 delta^2,gamma F_m)+V. The first channel has an instantaneous inverse and the second the new L1 scalar kernel. A weighted Volterra argument gives a unique prepared causal inverse at every fixed finite positive L, with continuous original lapse and scale reconstructed from continuous original forcing and smooth solutions for smooth prepared forcing. No numerical inverse norm or stability estimate is inferred.",
            "not_established":["Arbitrary spatial or independently varied initial-state/initial-jet response",
                "A numerical inverse norm, small response norm, quantum mode stability or causal cone bound",
                "A nonlinear neighborhood, global perturbed bounce or quantum response on a new background",
                "Canonical fixed-profile interactions, mixed loops, omitted operators, thresholds or Wilsonian cutoff",
                "Healthy full-candidate Minkowski vacuum, common-parent V/G/B, UV completion/exclusion or original P8 closure"],
            "verification_boundary":"Exact new time covariance, subtraction and finite-local Ward identities; literal tree-action reduction with both third-derivative cross channels; new dimensional fourth-contact and scalar normalization checks; written distributional matching and weighted Volterra proof with stronger local first-row reconstruction. Every parent is natively rebuilt. Finite analytic majorants are proved to exist but not numerically evaluated. Not proof-assistant formalization."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The prepared ordinary Proca inverse report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.87.PREPARED_ORDINARY_PROCA_INVERSE replay passed; quantum stability, matching and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
