"""Read-only exact massive first-sheet and causal-inverse certificate."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_vector_acoustic_response import verify as parent

from . import block, cut, kernel, proofs

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"massive-causal-inverse.json"
PARENT_SHA="8c1f2be8713fee678e977a07941510692cfd2f98a984b36db0e6cb180dbeb4c3"
sha,serialize=affine.sha,affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_vector_dispersion_inverse/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen acoustic covariance-response report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_71_fully_rebuilt":PARENT_SHA,"exact_massive_S6_69_block_not_replaced_by_leading_log":True}


def controls():
    calls=[]
    for value in (True,False,sp.true,sp.false,0,-1,1.0,sp.Float(1),"1",None,
                  sp.oo,-sp.oo,sp.nan,sp.zoo,sp.I,sp.Symbol("foreign"),[1]):
        calls.extend((lambda value=value:block.scale(value),
                      lambda value=value:kernel.moment_bound(value)))
    calls.extend(lambda value=value:serialize(value) for value in (1.0,sp.Float(1),sp.oo,sp.nan))
    rejected=0
    for call in calls:
        try:
            call()
        except (TypeError,ValueError):
            rejected+=1
    if rejected!=len(calls):
        raise ValueError("An unsupported massive inverse input was accepted")
    return {"rejected_inputs":rejected,"native_validation_before_symbolic_conversion":True,
            "mass_is_strictly_positive_no_massless_continuation":True,
            "positive_density_sign_not_reversed":True,
            "nonzero_instantaneous_matrix_not_deleted":True,
            "oscillatory_cut_integral_not_replaced_by_absolute_frequency_integral":True,
            "leading_log_spurious_root_not_used_as_an_exact_pole":True,
            "causal_distribution_inverse_not_a_generic_C1_endomorphism":True,
            "actual_curved_remainder_and_full_tree_plus_loop_inverse_not_claimed":True}


@cache
def build_report():
    prior=prior_checks()
    identities,gates=proofs.residuals(),proofs.checks()
    if not all(value is True for value in gates.values()):
        raise ValueError("A massive inverse audit gate failed")
    certified=affine.certify_residuals(identities)
    return {"schema":1,"claim":"P8-S6.72.MASSIVE_CAUSAL_INVERSE","date":"2026-09-08",
            "status":"EXACT_ISOLATED_FLAT_MASSIVE_BLOCK_CAUSAL_L1_INVERSE; ORIGINAL_P8_OPEN",
            "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
            "formulation":"FORMULATION.md",
            "written_proofs":["notes/first-sheet.md","notes/causal-kernel.md","notes/next-curved.md"],
            "exact_residuals":certified,"named_exact_check_count":len(identities),
            "checked_scalar_entries":sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in identities.values()),
            "proof_checks":gates,"first_sheet_block_and_endpoint_matrices":serialize(block.data()),
            "exact_cut_density_and_positive_Gram":serialize(cut.data()),
            "causal_kernel_endpoint_constants_and_representation":serialize(kernel.data()),
            "prepared_primitive_examples":serialize({"h_one":kernel.moment_bound(1),
                "h_interval_maximum":kernel.moment_bound(sp.Rational(125,64))}),
            "function_space_boundary":"Rinf times the identity plus a matrix L1 causal convolution is bounded on C0 and on bounded continuous half-line inputs. The exact positive moment also bounds its primitive. No numerical L1 norm is supplied; only the regular part has norm tending to zero on short intervals.",
            "variable_h_policy":"The exact matrix congruence defines S(t) B_1 S(t) with the inner multiplication inside the convolution. Its inverse is S(t)^-1 Rop_1 S(t)^-1. This is a defined preconditioner, not a derived curved-response identity.",
            "controls":controls(),
            "verdict":"The exact isolated flat massive dispersion block is zero-free on its first sheet and both cut banks. The inverse has a positive cut density, a fully determined nonzero rank-one instantaneous matrix and an integrable regular retarded kernel. This advances from real-axis matrix invertibility to an actual causal C0 inverse of that block only.",
            "not_established":["The actual curved selected-state no-loss stress remainder, full tree-plus-loop causal inverse, spatial or independent-state response",
                "Quantum stability/cones, interactions, other loops, cutoff, finite Wilson matching, common-parent V/G/B or original P8 closure"],
            "verification_boundary":"Exact symbolic matrix, sign, cut, endpoint and moment anchors with written first-sheet, keyhole contour, Abel-limit and time-endpoint proofs. Not proof-assistant formalization."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The massive causal-inverse report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.72.MASSIVE_CAUSAL_INVERSE replay passed; actual curved inverse and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
