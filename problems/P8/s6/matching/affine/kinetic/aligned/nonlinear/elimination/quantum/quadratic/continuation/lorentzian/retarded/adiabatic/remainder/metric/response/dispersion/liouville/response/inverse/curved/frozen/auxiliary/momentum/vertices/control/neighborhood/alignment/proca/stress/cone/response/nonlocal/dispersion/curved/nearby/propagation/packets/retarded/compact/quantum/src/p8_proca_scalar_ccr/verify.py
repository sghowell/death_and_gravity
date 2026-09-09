"""Read-only actual nearby reduced quadratic scalar CCR and Kubo-response certificate."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as original
from p8_proca_compact_response import verify as parent

from . import audit, canonical, growth, response, state

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"nearby-quadratic-scalar-CCR.json"
PARENT_SHA="b0980384d783fdb9b3cdbb4d8ee292a37011344ed01b313b0b1faaad1633e3a3"
sha,serialize=original.sha,original.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_proca_scalar_ccr/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen quantitative compact classical-response report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_92_fully_rebuilt":PARENT_SHA,
            "same_actual_nearby_background_quadratic_action_observable_and_compact_probes":True,
            "overall_positive_action_scale_and_hbar_retained":True,
            "new_scalar_Gaussian_state_not_a_transfer_of_old_Proca_quantum_profile":True}


@cache
def build_report():
    prior=prior_checks()
    identities,gates=audit.residuals(),audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("An actual nearby quadratic scalar CCR gate failed")
    exact=original.certify_residuals(identities)
    return {"schema":1,"claim":"P8-S6.93.NEARBY_QUADRATIC_SCALAR_CCR","date":"2026-09-09",
            "status":"ACTUAL_NEARBY_REDUCED_QUADRATIC_SCALAR_CCR_HAS_POSITIVE_TEMPERED_GAUSSIAN_LINEAR_FIELDS_AND_STATE_INDEPENDENT_COMPACT_MATTER_SPACELIKE_COMMUTATOR; HADAMARD_STRESS_INTERACTING_MATCHING_UV_AND_ORIGINAL_P8_OPEN",
            "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
            "formulation":"FORMULATION.md",
            "written_proofs":["notes/canonical.md","notes/growth.md","notes/state.md","notes/response.md","notes/scope.md"],
            "exact_residuals":exact,"named_exact_check_count":len(identities),
            "checked_scalar_entries":sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in identities.values()),
            "proof_checks":gates,
            "literal_density_Hamiltonian_symplectic_map_and_source_normalization":serialize(canonical.data()),
            "actual_all_momentum_transfer_growth_and_Schwartz_test_space":serialize(growth.data()),
            "explicit_positive_Gaussian_Cauchy_covariance_and_field_state":serialize(state.data()),
            "state_independent_Kubo_response_and_compact_commutator_lower":serialize(response.data()),
            "controls":audit.controls(),
            "verdict":"The specified reduced quadratic scalar Hamiltonian on the actual nearby local classical bounce has canonical density evolution preserving the exact CCR. Complete high- and low-frequency transfer bounds and a written all-orders parameter-ODE induction make evolved Schwartz and compact spacetime smearings well-defined. An explicit initial positive Gram covariance constructs a tempered Gaussian linear-field state. The normalized linear Kubo response equals the classical retarded response in every state of the same CCR algebra. The unchanged S6.92 compact probes yield a central matter-spacelike commutator with magnitude at least 3*hbar*D/(4*kappa)>0, keeping the action scale and source volume. This is not a Hadamard or semiclassical state, an interacting response, a cutoff-valid UV exclusion or a nonperturbative quantum-gravity locality theorem.",
            "not_established":["Hadamard singularities, finite state energy, local renormalized Wick polynomials, stress tensor or nearby semiclassical backreaction",
                "A global infinite-volume Fock-unitary implementation of time evolution",
                "A transfer of the old Proca state, background quantum profile, stress or physical scale dictionary",
                "Canonical interactions, omitted operators, loop or heavy-threshold errors, joint temporal/spatial EFT validity or actual parent high-band tail",
                "Finite Wilson matching, healthy full-candidate vacuum, nonperturbative local quantum-gravity observables or common-parent V/G/B",
                "A globally complete nearby bounce, universal UV completion/exclusion, realistic-field P8(a) or original P8 closure"],
            "verification_boundary":"Native exact canonical density Hamiltonian and CCR identities, complete endpoint map and transfer majorants, explicit positive covariance Gram factor and antisymmetric part, Kubo sign/source/action normalization and inherited compact response lower, with written all-orders Schwartz, Fock-state and linear-source proofs. Every parent is rebuilt natively. No scientific-library patch, finite-jet substitute for test-space continuity, transplanted Klein-Gordon theorem, Hadamard/stress assertion or interacting UV transfer. Not proof-assistant formalization."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The nearby quadratic scalar CCR report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.93.NEARBY_QUADRATIC_SCALAR_CCR replay passed; Hadamard stress, interacting matching, UV and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
