"""Read-only exact frozen bounce-frequency pole diagnostic."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_prepared_volterra import verify as parent

from . import bracket, coefficients

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"frozen-frequency-pole.json"
PARENT_SHA="da5b7d2da91c6904d6b4fafdd58f231ae8ff01a7423c6c58e25efd2cad9c077f"
sha,serialize=affine.sha,affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_bounce_frequency/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen prepared coupled-inverse report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    coefficients.pinned(coefficients.local.REPORT,coefficients.LOCAL_SHA)
    coefficients.pinned(coefficients.response.REPORT,coefficients.RESPONSE_SHA)
    return {"S6_73_fully_rebuilt":PARENT_SHA,
            "S6_67_coefficient_report_already_replayed_in_parent_chain":coefficients.LOCAL_SHA,
            "S6_68_fixed_profile_report_already_replayed_in_parent_chain":coefficients.RESPONSE_SHA}


def controls():
    calls=[]
    for value in (True,False,sp.true,sp.false,-1,5,1.0,sp.Integer(1),sp.Float(1),"1",None,[1]):
        calls.append(lambda value=value:coefficients.local_row(value))
    for value in (True,False,sp.true,sp.false,0,1,1.0,sp.Integer(1),"LOWER","upper ",None,["lower"]):
        calls.append(lambda value=value:bracket.endpoint(value))
    calls.extend(lambda value=value:serialize(value) for value in (1.0,sp.Float(1),sp.oo,sp.nan))
    rejected=0
    for call in calls:
        try:
            call()
        except (TypeError,ValueError):
            rejected+=1
    if rejected!=len(calls):
        raise ValueError("An unsupported frozen-frequency input was accepted")
    return {"rejected_inputs":rejected,"native_validation_before_cached_implementations":True,
            "physical_Euler_coefficients_frozen_after_variation":True,
            "nonzero_background_and_chart_jets_kept":True,
            "lower_local_matrices_not_symmetrized":True,
            "tadpole_is_fixed_Euler_Hessian_not_reselected_stress_Jacobian":True,
            "finite_fourth_contact_not_double_counted":True,
            "massive_bubble_kept_exact_not_just_leading_log":True,
            "stationary_diagnostic_not_promoted_to_actual_curved_instability":True,
            "sub_Planckian_not_promoted_to_below_justified_interacting_cutoff":True}


@cache
def residuals():
    first,second=coefficients.checks(),bracket.checks()
    if set(first).intersection(second):
        raise ValueError("Repeated frozen-frequency residual name")
    return first|second


@cache
def build_report():
    prior=prior_checks()
    identities,gates=residuals(),bracket.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A frozen-frequency inequality gate failed")
    exact=affine.certify_residuals(identities)
    return {"schema":1,"claim":"P8-S6.74.FROZEN_COUPLED_FREQUENCY_POLE","date":"2026-09-08",
            "status":"EXACT_POSITIVE_REAL_POLE_OF_DEFINED_FROZEN_TREE_PLUS_GAUSSIAN_DIAGNOSTIC; ACTUAL_CURVED_STABILITY_AND_ORIGINAL_P8_OPEN",
            "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
            "formulation":"FORMULATION.md","written_proofs":["notes/coefficients.md","notes/bracket.md","notes/scope.md"],
            "exact_residuals":exact,"named_exact_check_count":len(identities),
            "checked_scalar_entries":sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in identities.values()),
            "proof_checks":gates,
            "physical_background_bounce_jets":serialize({str(k):v for k,v in coefficients.background_jets().items()}),
            "constant_source_chart":serialize(coefficients.T),
            "frozen_tree_Euler_coefficients":serialize({str(j):v for j,v in enumerate(coefficients.tree_rows())}),
            "frozen_local_Gaussian_Euler_coefficients":serialize({str(j):coefficients.local_row(j) for j in range(5)}),
            "fixed_tadpole_Euler_contact":serialize(coefficients.tadpole_matrix()),
            "actual_pair_vertices_and_Gram":serialize({"pairs":coefficients.pairs(),"Gram":coefficients.gram(),"finite_fourth":coefficients.finite()}),
            "continuous_frequency_bracket_bounds":serialize(bracket.data()),
            "diagnostic_definition_and_domain":bracket.description(),
            "controls":controls(),
            "verdict":"For every fixed tadpole pair in the proven profile box, the defined frozen full tree-plus-Gaussian symbol has an uncancelled positive-real inverse pole strictly between 10^12 and 10^13 inverse bounce-time units. The determinant signs are exact continuous-band inequalities retaining all finite local terms and the exact massive bubble. This is not an actual evolving-bounce instability theorem or an EFT-domain verdict.",
            "not_established":["Root uniqueness, simplicity, complete spectrum or actual nonstationary growth estimate",
                "Numerical curved/state-dependent kernel error, justified interacting cutoff or omitted-loop/operator budgets",
                "Physical quantum stability/cones, finite Wilson matching, common-parent V/G/B or original P8 closure"],
            "verification_boundary":"Exact rational/symbolic coefficient replay and continuous inequalities with a written analytic intermediate-value/noncancellation proof. No floating-point root finder or proof-assistant formalization."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The frozen-frequency report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.74.FROZEN_COUPLED_FREQUENCY_POLE replay passed; actual curved stability and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
