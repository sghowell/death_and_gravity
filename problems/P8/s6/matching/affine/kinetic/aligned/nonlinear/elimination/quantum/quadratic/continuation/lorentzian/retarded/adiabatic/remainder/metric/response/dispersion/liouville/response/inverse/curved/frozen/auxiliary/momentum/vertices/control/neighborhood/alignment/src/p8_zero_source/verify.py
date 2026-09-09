"""Read-only separately named zero-source action and finite-change certificate."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as original
from p8_offclock_scalar import verify as parent

from . import affine, audit, bounds, coefficients, model, neighborhood, response, tree

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"zero-source-candidate.json"
PARENT_SHA="655ccb27216ec50271ca2da704a4d424f52380b23d37ad93590ffc13603a93fc"
sha,serialize=original.sha,original.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_zero_source/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen off-clock coupled-symbol report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_79_fully_rebuilt":PARENT_SHA,
            "old_action_negative_diagnostic_and_on_clock_tree_domain_not_overwritten":True}


@cache
def build_report():
    prior=prior_checks()
    identities,gates=audit.residuals(),audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A new zero-source action gate failed")
    exact=original.certify_residuals(identities)
    m,a=model.data(),affine.data()
    return {"schema":1,"claim":"P8-S6.80.ZERO_SOURCE_CANDIDATE","date":"2026-09-08",
            "status":"NEW_LOCAL_AFFINE_SOURCE_CONTACT_CANDIDATE_WITH_SAME_QUADRATIC_THEORY_AND_FINITE_HARD_TREE_CHANGE_BUDGET; NEIGHBORHOOD_HEALTH_MATCHING_AND_ORIGINAL_P8_OPEN",
            "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
            "formulation":"FORMULATION.md",
            "written_proofs":["notes/action-affine.md","notes/auxiliary.md","notes/finite-budget.md",
                              "notes/response-neighborhood.md","notes/scope.md"],
            "exact_residuals":exact,"named_exact_check_count":len(identities),
            "checked_scalar_entries":sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in identities.values()),
            "proof_checks":gates,
            "literal_new_action_and_Hamiltonian":serialize({
                "new_Hamiltonian":m["new_Hamiltonian"],"new_minus_old_Hamiltonian":m["new_minus_old_Hamiltonian"],
                "new_trace_reconstruction":m["new_trace_reconstruction"],
                "new_temporal_reconstruction":m["new_temporal_reconstruction"],
                "actual_original_source_jets":model.source_jets(),
                "literal_added_density":"N*sqrt(h_hat)*U*(T*S_normal-S_normal^2/2)/gamma_t"}),
            "all_connection_affine_lift":serialize({
                "added_connection_source":a["added_connection_source"],
                "unchanged_retained_mass_matrix":a["unchanged_retained_mass_matrix"],
                "normal_source_vector":a["normal_source_vector"],
                "zero_full_connection_Hessian_change":True,
                "all64_stationary_Euler_and_projective_complement_checks_in_residual_inventory":True}),
            "actual_new_five_lapse_jets_for_fourteen_inventory_rows":serialize(coefficients.jets()),
            "continuous_new_auxiliary_domain":serialize(bounds.auxiliary()),
            "universal_new_coefficient_bounds_and_positive_scaling":serialize({
                "coefficient_comparison":bounds.coefficients(),"degree_proof":bounds.scaling()}),
            "finite_same_scale_hard_tree_change":serialize(tree.data()),
            "fixed_light_Gaussian_and_changed_source_functional":serialize(response.data()),
            "old_offclock_mechanism_comparison":serialize(neighborhood.data()),
            "canonical_comparison_convention":"A literal new action, not a field redefinition. The two specified nonlinear canonical phase charts are compared by their stated identity, with identical tangent free variables and quadratic Hamiltonian. Old and new nonlinear momenta are not identified as the same Legendre map.",
            "controls":audit.controls(),
            "verdict":"The new local normal source/contact change lifts to the original affine system with unchanged connection Hessian and projective kernel. It preserves the old clock and complete quadratic theory, with unchanged original vector covariance and fixed-light Gaussian determinant. Actual new coefficient jets and continuous auxiliary bounds give a finite H3/H4 plus H3-squared model-change budget at the separately named M*tau=10^400. The particular old off-clock scalar-longitudinal q^4 mechanism is removed. The off-clock mass ratio, same-scale quadratic quantum diagnostic and unbounded matching requirements are not cured or hidden.",
            "not_established":["Full physical neighborhood cone, all-frequency or actual time-dependent stability",
                "An interacting Wilsonian cutoff, omitted operators, loops or heavy-threshold error bounds",
                "Transfer of the old L=10^24 numerical pole bracket to the new L=10^400 example",
                "A new off-clock Hadamard theorem or a bound on mixed vector-light/graviton loops",
                "Finite Wilson matching, common-parent V/G/B, UV completion/exclusion or original P8 closure"],
            "verification_boundary":"Exact local affine/ADM/coefficient identities, rational enclosures and written positive-majorant, Gaussian-transfer and mechanism-comparison proofs. The prior full chain is natively rebuilt. Not proof-assistant formalization."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The zero-source candidate report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.80.ZERO_SOURCE_CANDIDATE replay passed; cutoff, matching and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
