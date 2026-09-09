"""Read-only constant-mass physical Proca candidate and quantum comparison."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as original
from p8_zero_source import verify as parent

from . import affine, audit, bounds, model, propagation, quantum

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"constant-physical-Proca.json"
PARENT_SHA="f33dce0fcfc685ebf174d475d458948dfdc81ac2c52179ad494cd05b3ebf1d65"
sha,serialize=original.sha,original.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_constant_proca/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen zero-source candidate report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_80_fully_rebuilt":PARENT_SHA,
            "previous_action_and_quantum_nontransfer_boundary_preserved":True}


@cache
def build_report():
    prior=prior_checks()
    identities,gates=audit.residuals(),audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A constant-mass Proca candidate gate failed")
    exact=original.certify_residuals(identities)
    a=affine.data()
    return {"schema":1,"claim":"P8-S6.81.CONSTANT_PHYSICAL_PROCA","date":"2026-09-08",
            "status":"NEW_CLASSICAL_AFFINE_CONSTANT_MASS_PROCA_CANDIDATE_WITH_PHYSICAL_VECTOR_CONE_AND_FINITE_TREE_CHANGE; QUANTUM_TADPOLE_RESPONSE_MATCHING_AND_ORIGINAL_P8_OPEN",
            "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
            "formulation":"FORMULATION.md",
            "written_proofs":["notes/affine.md","notes/classical.md","notes/propagation.md",
                              "notes/quantum.md","notes/scope.md"],
            "exact_residuals":exact,"named_exact_check_count":len(identities),
            "checked_scalar_entries":sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in identities.values()),
            "proof_checks":gates,
            "literal_affine_mass_source_and_contact_update":serialize({
                "old_retained_mass":a["old_retained_mass"],"new_retained_mass":a["new_retained_mass"],
                "retained_mass_change":a["retained_mass_change"],
                "original_one_form_shift":a["original_one_form_shift"],
                "literal_added_affine_density":a["literal_added_affine_density"],
                "added_full_connection_source":a["added_full_connection_source"],
                "quotient_determinant_ratio":a["quotient_determinant_ratio"],
                "connection_Hessian_changes":True,
                "full64_Euler_inverse_lift_projective_and_complement_checks_in_residuals":True}),
            "actual_new_Hamiltonian":serialize(model.data()),
            "actual_new_lapse_jet_inventory":serialize(model.jets()),
            "continuous_new_auxiliary_and_coefficient_bounds":serialize(bounds.auxiliary()),
            "same_scale_finite_hard_tree_change":serialize(bounds.tree()),
            "constrained_physical_Proca_symbol":serialize(propagation.symbol()),
            "actual_nearby_vector_canonical_blocks":serialize(propagation.neighborhood()),
            "changed_quantum_readouts_and_operator_jets":serialize(quantum.matrices()),
            "new_constant_mass_local_quantum_coefficients_and_boundary":serialize(quantum.boundary()),
            "controls":audit.controls(),
            "verdict":"The literal new affine mass/source/contact update gives an invertible local quotient and a source-free constant positive physical Proca mass. The original classical clock and full quadratic seven-mode theory agree, with continuous auxiliary and finite hard-tree change bounds on the same L=10^400 domain. The constrained vector characteristic cone is the physical metric cone on admitted smooth globally hyperbolic backgrounds; no full scalar-neighborhood or nonlinear result is inferred. Actual quantum energy readouts and first/second lapse operator jets change, although original clock modes and covariance agree. The old fixed tadpole cancellation and quantum response are not transferred.",
            "not_established":["New finite stress/tadpole derivative bounds and exact quantum-corrected clock cancellation",
                "The new Gaussian light response, its spectrum or all-frequency quantum stability",
                "Full coupled scalar neighborhood cone, nonlinear spatial well-posedness or actual global stability",
                "Interacting Wilsonian cutoff, omitted operators, mixed loops or heavy-threshold errors",
                "Finite Wilson matching, common-parent V/G/B, UV completion/exclusion or original P8 closure"],
            "verification_boundary":"Exact affine/ADM/coefficient/symbol/contact identities, rational enclosures and written proofs. Vector finite propagation uses the stated constrained-Proca Green-hyperbolic theorem with its hypotheses, not an unconstrained four-mode claim. The prior full chain is natively rebuilt. Not proof-assistant formalization."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The constant-mass Proca report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.81.CONSTANT_PHYSICAL_PROCA replay passed; quantum matching and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
