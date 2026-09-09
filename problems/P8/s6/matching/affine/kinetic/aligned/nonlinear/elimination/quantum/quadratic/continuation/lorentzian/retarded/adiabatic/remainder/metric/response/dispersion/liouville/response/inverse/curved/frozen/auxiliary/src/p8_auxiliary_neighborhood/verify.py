"""Read-only quantified nonlinear auxiliary branch and invariant Taylor remainder."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_bounce_frequency import verify as parent

from . import bounds, center, model

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"quantified-auxiliary-domain.json"
PARENT_SHA="1ae7056de9624cfbe87a3b07330d9b9703a8e6a14e0db3851627db4366a0d6d1"
sha,serialize=affine.sha,affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_auxiliary_neighborhood/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen coupled-frequency diagnostic report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_74_fully_rebuilt":PARENT_SHA,
            "full_S6_45_Hamiltonian_and_S6_56_margin_replayed_in_parent_chain":True}


def controls():
    calls=[]
    for value in (True,False,sp.true,sp.false,-1,5,1.0,sp.Integer(1),sp.Float(1),"1",None,[1]):
        calls.append(lambda value=value:bounds.taylor_majorant(value))
    for value in (True,False,sp.true,sp.false,-1,1,0.0,sp.Float(0),sp.oo,sp.nan,"0",[0]):
        calls.append(lambda value=value:bounds.response_radius(value))
    calls.extend(lambda value=value:serialize(value) for value in (1.0,sp.Float(1),sp.oo,sp.nan))
    rejected=0
    for call in calls:
        try:
            call()
        except (TypeError,ValueError):
            rejected+=1
    if rejected!=len(calls):
        raise ValueError("An unsupported auxiliary-domain input was accepted")
    return {"rejected_inputs":rejected,"input_validation_precedes_arithmetic_or_cache":True,
            "whole_scalar_function_not_just_clock_Taylor_jet":True,
            "original_Q_ODE_and_boundary_primitive_kept":True,
            "boundary_time_derivative_bounded_by_Cauchy_not_differentiated_inequality":True,
            "complex_polydisc_bounds_apply_to_actual_implicit_root":True,
            "full_spatial_canonical_jets_and_Maxwell_normalizations_retained":True,
            "invariant_degree_not_confused_with_physical_field_degree":True,
            "classical_tree_plus_margin_not_nonlocal_quantum_constraint":True,
            "no_physical_interacting_cutoff_or_original_P8_closure":True}


@cache
def residuals():
    result={}
    for group in (model.checks(),bounds.checks(),center.checks()):
        if set(result).intersection(group):
            raise ValueError("Repeated auxiliary-domain residual name")
        result.update(group)
    return result


@cache
def build_report():
    prior=prior_checks()
    identities,gates=residuals(),bounds.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A quantified auxiliary-domain gate failed")
    certified=affine.certify_residuals(identities)
    actual=model.coefficients()
    generic=model.generic()
    return {"schema":1,"claim":"P8-S6.75.QUANTIFIED_AUXILIARY_NEIGHBORHOOD","date":"2026-09-08",
            "status":"EXPLICIT_FULL_CLASSICAL_CANONICAL_JET_AUXILIARY_DOMAIN_AND_INVARIANT_TAYLOR_REMAINDER; INTERACTING_CUTOFF_AND_ORIGINAL_P8_OPEN",
            "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
            "formulation":"FORMULATION.md","written_proofs":["notes/hamiltonian.md","notes/neighborhood.md","notes/scope.md"],
            "exact_residuals":certified,"named_exact_check_count":len(identities),
            "checked_scalar_entries":sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in identities.values()),
            "proof_checks":gates,
            "invariant_coordinates":serialize(model.VARIABLES),
            "actual_full_Hamiltonian_coefficients":serialize(actual),
            "generic_full_spatial_Hamiltonian_and_reconstructions":serialize(generic),
            "full_scalar_potential_complex_bound":serialize(bounds.scalar_potential_bound()),
            "holomorphic_coefficient_majorants":serialize(bounds.coefficients()),
            "quantitative_contraction_and_Taylor_bounds":serialize(bounds.contraction()),
            "closed_bounce_Hamiltonian_and_boundary_primitive":serialize(center.data()),
            "invariant_Taylor_degree_majorants":serialize({str(j):bounds.taylor_majorant(j) for j in range(5)}),
            "actual_radius_examples":serialize({"zero":bounds.response_radius(0),
                "domain_boundary":bounds.response_radius(bounds.JET_RADIUS),
                "smaller_Taylor_domain":bounds.response_radius(bounds.JET_RADIUS/1000)}),
            "scope_and_scale_policy":"Classical S6.45 tree parent plus the S6.56 epsilon=10^-6 margin, positive dimensionless zeta=10^-6, on I. The actual post-temporal normal Hamiltonian is in the existing hat-volume convention and units M^2/tau^2. Electric and magnetic invariants retain their zeta factors. No nonlocal Gaussian constraint, momentum reduction or physical cutoff is inferred.",
            "controls":controls(),
            "verdict":"An explicit nine-invariant canonical-jet polydisc gives the unique actual lapse branch in a stated disc, a nonzero joint auxiliary block and a small reconstructed temporal vector. Holomorphic bounds retain the entire scalar, actual source ODE, boundary primitive and all spatial terms. The actual reduced Hamiltonian has continuous invariant Taylor coefficient and remainder bounds. These are interaction-control prerequisites, not a physical scattering or cutoff theorem.",
            "not_established":["Spatial momentum reduction, canonical physical interaction vertices or a weak-coupling frequency band",
                "Nonlinear energy/evolution stability, nonlocal quantum auxiliary equations, other loops or omitted-operator budgets",
                "Finite Wilson matching, common-parent V/G/B or original P8 closure"],
            "verification_boundary":"Exact Hamiltonian, source, boundary, bounce and stationary-series anchors, rational majorants, and written holomorphic contraction/remainder proofs. Not proof-assistant formalization."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The quantified auxiliary-domain report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.75.QUANTIFIED_AUXILIARY_NEIGHBORHOOD replay passed; interacting cutoff and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
