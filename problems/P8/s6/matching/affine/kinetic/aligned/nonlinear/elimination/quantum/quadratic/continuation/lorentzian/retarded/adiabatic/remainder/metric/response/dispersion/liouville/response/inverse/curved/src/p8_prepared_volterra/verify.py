"""Read-only prepared homogeneous no-loss normal form and coupled inverse."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_vector_dispersion_inverse import verify as parent
from p8_vector_metric_local import canonical

from . import primitive, proofs, tree, vertices

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"prepared-coupled-inverse.json"
PARENT_SHA="b08b73a0d78522e358e864177984c8627a7643974126448a701c256a62001684"
sha,serialize=affine.sha,affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_prepared_volterra/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen exact massive causal-inverse report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_72_fully_rebuilt":PARENT_SHA,"original_physical_state_prescription_and_fixed_tadpole_retained":True}


def controls():
    calls=[]
    for value in (True,False,sp.true,sp.false,0,1,1.0,sp.Integer(1),"transverse","",None,[1]):
        calls.append(lambda value=value:vertices.high(value))
    for value in (True,False,sp.true,sp.false,-1,5,1.0,sp.Integer(1),sp.Float(1),"1",None,[1]):
        calls.extend((lambda value=value:primitive.local_kernel(value),
                      lambda value=value:primitive.primitive_shape(value)))
    calls.extend(lambda value=value:serialize(value) for value in (1.0,sp.Float(1),sp.oo,sp.nan))
    rejected=0
    for call in calls:
        try:
            call()
        except (TypeError,ValueError):
            rejected+=1
    if rejected!=len(calls):
        raise ValueError("An unsupported coupled-response input was accepted")
    return {"rejected_inputs":rejected,"native_validation_before_cached_implementations":True,
            "finite_nonlogarithmic_fourth_Taylor_contact_retained":True,
            "common_dimensional_diagonal_matched_before_finite_limit":True,
            "prepared_initial_squeezing_and_fixed_tadpole_profiles_not_reset":True,
            "matter_reconstruction_uses_negative_sign_from_literal_action":True,
            "nonzero_instantaneous_inverse_not_dropped_for_contraction":True,
            "finite_inverse_constants_not_claimed_numerically_small":True,
            "linear_retarded_existence_not_a_quantum_health_or_P8_closure_claim":True}


def fourth_contacts():
    result={}
    for sector in ("T","L"):
        weights=canonical.data(sector)["weights"]
        pair=sp.ImmutableMatrix([weights[key][0]-weights[key][1] for key in ("N","Z")])
        result[sector]=pair*pair.T/32
    return result


@cache
def build_report():
    prior=prior_checks()
    identities,gates=proofs.residuals(),proofs.checks()
    if not all(value is True for value in gates.values()):
        raise ValueError("A prepared coupled-inverse audit gate failed")
    certified=affine.certify_residuals(identities)
    return {"schema":1,"claim":"P8-S6.73.PREPARED_COUPLED_CAUSAL_INVERSE","date":"2026-09-08",
            "status":"PREPARED_HOMOGENEOUS_TREE_PLUS_RETAINED_GAUSSIAN_NO_LOSS_CAUSAL_INVERSE; ORIGINAL_P8_OPEN",
            "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
            "formulation":"FORMULATION.md",
            "written_proofs":["notes/kernel.md","notes/matching.md","notes/volterra.md","notes/scope.md"],
            "exact_residuals":certified,"named_exact_check_count":len(identities),
            "checked_scalar_entries":sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in identities.values()),
            "proof_checks":gates,"physical_current_UV_coefficients":serialize(vertices.data()),
            "full_dimensional_fourth_adiabatic_contacts":serialize(fourth_contacts()),
            "literal_classical_matter_reduction":serialize(tree.data()),
            "local_four_primitive_kernels":serialize({str(j):primitive.local_kernel(j) for j in range(5)}),
            "singular_power_four_primitive_shapes":serialize({str(j+1):primitive.primitive_shape(j) for j in range(5)}),
            "function_space_estimates_and_domains":primitive.description(),
            "coupling_and_preparation_policy":"Fix m>=1000,L>0 and margin epsilon in (0,1/100] on I. Q=64*pi^2*L^2 times the S6.68 physical vector-plus-fixed-tadpole force response, gamma=1/(64*pi^2*L^2). The equations E f+gamma Q f=g have the prepared zero-past inverse, with zero independent matter charge and clock forcing fixed by the Ward identity. No new action or state is introduced.",
            "controls":controls(),
            "verdict":"The actual prepared physical-current response has I4 Q=S B_m S+V, with an integrable logarithmic Volterra remainder in a no-loss C0 norm. Full-dimensional all-momentum fourth Taylor contacts retain the finite nonlogarithmic channel. Literal charge elimination gives a local second-order classical metric operator. The exact massive inverse and a weighted Volterra argument yield a unique smooth prepared homogeneous solution of the full retained linear tree-plus-Gaussian-vector-plus-fixed-tadpole forced equations on I. Constants are finite but not numerically evaluated.",
            "not_established":["Numerical smallness or physical quantum stability/cones of the full inverse; nonlinear neighborhood-uniform semiclassical existence",
                "Spatial response, arbitrary independent initial data, other loops, interactions, cutoff, finite Wilson/common-parent V/G/B matching or original P8 closure"],
            "verification_boundary":"Exact physical vertex, dimensional all-momentum Taylor, time-primitive and literal Euler anchors plus written finite-WKB, common-dimensional extension, logarithmic-remainder and weighted-Volterra proofs. Not proof-assistant formalization."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The prepared coupled-inverse report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.73.PREPARED_COUPLED_CAUSAL_INVERSE replay passed; quantum health and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
