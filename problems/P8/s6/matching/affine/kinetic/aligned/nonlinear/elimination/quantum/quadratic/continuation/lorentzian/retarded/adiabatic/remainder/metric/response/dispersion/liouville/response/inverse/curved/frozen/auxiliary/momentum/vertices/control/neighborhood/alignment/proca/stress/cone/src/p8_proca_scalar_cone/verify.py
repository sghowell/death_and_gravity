"""Read-only actual central scalar cone and finite-amplitude boundary."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as original
from p8_proca_stress import verify as parent

from . import audit, background, domain, principal, spatial

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"central-scalar-cone.json"
PARENT_SHA="4e3c6103906b8ddb344806b7fe9309493f3c8dfda28629e70912872be65c09f8"
sha,serialize=original.sha,original.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_proca_scalar_cone/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen ordinary Proca stress/profile report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_82_fully_rebuilt":PARENT_SHA,
            "classical_S6_81_action_is_distinguished_from_S6_82_fixed_quantum_profile":True}


@cache
def build_report():
    prior=prior_checks()
    identities,gates=audit.residuals(),audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("An actual central scalar cone gate failed")
    exact=original.certify_residuals(identities)
    return {"schema":1,"claim":"P8-S6.83.CENTRAL_SCALAR_CONE","date":"2026-09-08",
            "status":"ACTUAL_EULER_FIRST_CLASSICAL_CENTRAL_SCALAR_CONE_AND_FINITE_AMPLITUDE_BOUNDARY; FULL_QUANTUM_NEIGHBORHOOD_MATCHING_AND_ORIGINAL_P8_OPEN",
            "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
            "formulation":"FORMULATION.md",
            "written_proofs":["notes/background.md","notes/spatial.md","notes/cones.md",
                              "notes/boundary.md","notes/scope.md"],
            "exact_residuals":exact,"named_exact_check_count":len(identities),
            "checked_scalar_entries":sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in identities.values()),
            "proof_checks":gates,
            "actual_homogeneous_background_time_jets":serialize(background.data()),
            "full_spatial_two_leg_time_jet_reduction":serialize({
                left+"__"+right:spatial.pair(left,right)
                for a,left in enumerate(spatial.CHANNELS) for right in spatial.CHANNELS[a:]}),
            "actual_canonical_and_Euler_first_blocks":serialize(principal.blocks()),
            "physical_principal_matrices_and_characteristic_factorization":serialize(principal.formula()),
            "continuous_central_positive_subluminal_domain":serialize(domain.domain()),
            "actual_finite_amplitude_luminality_boundary_and_superluminal_control":serialize(domain.boundary()),
            "controls":audit.controls(),
            "verdict":"For the source-free constant-physical-mass Proca classical candidate, all scalar phase pairs and their first background time jets are reduced before taking the principal limit. On the central homogeneous family |N-1|<=10^-8, both scalar kinetic and gradient forms are positive, the clock is strictly subluminal and matter is exactly luminal in the physical metric. The positive tensor and ordinary Proca principal cones remain luminal. The same classical family crosses clock luminality at a unique simple root between N=1+2*10^-7 and N=1+2.1*10^-7; the exact datum N=1+10^-6 is ghost/gradient positive but superluminal. No old off-clock vector-mixing argument or old interaction domain is transferred.",
            "not_established":["Uniform nonlinear or inhomogeneous neighborhood well-posedness, or a global nearby bounce",
                "Scalar cones after adding the fixed quantum profile and full Gaussian response",
                "A finite-frequency superluminal advance inside an independently controlled interacting EFT band",
                "Old nine-invariant interaction-polydisc or finite hard-tree bounds on this broader central interval",
                "Mixed loops, omitted operators, V/G/B matching, UV completion/exclusion or original P8 closure"],
            "verification_boundary":"Exact canonical background equations, full two-leg spatial reduction, Euler-first principal matrices, rational factorization and continuous interval enclosures with explicit negative controls. Parent reports are natively rebuilt. A classical central-slice family only; not proof-assistant formalization or a full quantum solution."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The central scalar cone report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.83.CENTRAL_SCALAR_CONE replay passed; full quantum neighborhood, matching and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()

