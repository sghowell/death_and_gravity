"""Read-only actual nearby local classical ordinary-Proca bounce solutions."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as original
from p8_proca_prepared_inverse import verify as parent

from . import audit, existence, jets, taylor

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"local-nearby-Proca-bounces.json"
PARENT_SHA="25eabbd3c7847a746ad7d9b015f27cf63cb3e23aa6d6acfb9b8697fc720bf8fc"
sha,serialize=original.sha,original.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_proca_nearby_bounce/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen prepared ordinary Proca coupled inverse report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_87_fully_rebuilt":PARENT_SHA,
            "new_classical_backgrounds_do_not_inherit_old_quantum_state_profile_or_inverse":True}


@cache
def build_report():
    prior=prior_checks()
    identities,gates=audit.residuals(),audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("An actual nearby classical bounce solution gate failed")
    exact=original.certify_residuals(identities)
    return {"schema":1,"claim":"P8-S6.88.LOCAL_NEARBY_PROCA_BOUNCES","date":"2026-09-09",
            "status":"ACTUAL_REAL_ANALYTIC_NEARBY_CLASSICAL_HOMOGENEOUS_BOUNCES_ON_UNIFORM_CLOCK_INTERVAL_ABS_U_LE_10_MINUS_7; TIME_DEPENDENT_CONES_QUANTUM_MATCHING_AND_ORIGINAL_P8_OPEN",
            "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
            "formulation":"FORMULATION.md",
            "written_proofs":["notes/acceleration.md","notes/complex-domain.md","notes/evolution.md",
                              "notes/enlarged-domain.md","notes/scope.md","notes/sources.md"],
            "exact_residuals":exact,"named_exact_check_count":len(identities),
            "checked_scalar_entries":sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in identities.values()),
            "proof_checks":gates,
            "actual_central_constraint_and_physical_acceleration":serialize(jets.data()),
            "continuous_central_family_enclosures":serialize(jets.boxes()),
            "exact_representative_central_data":serialize({str(value):existence.datum(value) for value in
                (sp.Integer(0),sp.Rational(1,10**8),sp.Rational(1,10**6))}),
            "fresh_joint_complex_coefficient_domain":serialize(existence.coefficient_domain()),
            "coarse_lapse_root":serialize(existence.lapse_root()),
            "coarse_evolution":serialize(existence.evolution()),
            "enlarged_literal_central_jet":serialize(taylor.literal_jet()),
            "enlarged_exact_jet_coefficients_and_enclosures":serialize(taylor.coefficients()),
            "enlarged_lapse_and_evolution_bounds":serialize(taylor.bounds()),
            "enlarged_time_continuations":serialize({str(value):taylor.evolution(value) for value in
                (taylor.TIME_WINDOW,taylor.TIME_WINDOW/2)}),
            "controls":audit.controls(),
            "verdict":"For every real central lapse shift in [0,10^-6], the S6.81 source-free classical action has a constraint-compatible real analytic homogeneous solution on |u|<=10^-7 with its fixed nonzero free-matter charge. A fresh complex-time Hamiltonian bound, exact central jets and radial Cauchy tails give a lapse-root contraction below 1/8 and a two-variable Picard contraction 2/5. Exact reconstruction and the Noether identity recover the original equations. The physical metric, not the hat scale alone, has H(0)=0, proper Hubble acceleration above 3.9999 and a(0)>=1. The third time derivative of the original primitive is retained in the lapse acceleration. These actual local bounces realize the S6.83 central cone data, including its superluminal control; no resolved propagation or UV verdict follows.",
            "not_established":["Global complete nearby bounces, all-time lower scale bounds or nonlinear inhomogeneous stability",
                "A quantitative time-dependent scalar cone, wavepacket advance or admitted propagation band",
                "The S6.82 profile-corrected quantum action, new quantum state, cancellation or response/inverse on these backgrounds",
                "Canonical interactions, mixed loops, omitted operators, thresholds or an interacting Wilsonian cutoff on the new domain",
                "Finite Wilson matching, common-parent V/G/B, UV completion/exclusion or original P8 closure"],
            "verification_boundary":"Exact fixed-phase Hamiltonian and primitive derivatives, independent canonical reduction, continuous central rational enclosures, actual finite Taylor jets and explicit analytic tail bounds, with written complex contraction, Picard, parity, reconstruction and Noether proofs. Every parent is natively rebuilt. The executable rational datum interface is a finite interface to a written theorem on the whole real central interval. Not proof-assistant formalization."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The local nearby Proca bounce report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.88.LOCAL_NEARBY_PROCA_BOUNCES replay passed; time-dependent cones, matching and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
