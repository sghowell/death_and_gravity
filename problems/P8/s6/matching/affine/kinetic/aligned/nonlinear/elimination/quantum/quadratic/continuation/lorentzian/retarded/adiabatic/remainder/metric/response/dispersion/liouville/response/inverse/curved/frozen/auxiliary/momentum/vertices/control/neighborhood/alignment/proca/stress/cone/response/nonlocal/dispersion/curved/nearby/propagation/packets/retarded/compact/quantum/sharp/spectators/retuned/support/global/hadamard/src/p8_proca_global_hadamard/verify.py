"""Read-only actual global two-cone Hadamard scalar-state certificate."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as certificate
from p8_proca_global_support import verify as parent

from . import (
    audit,
    canonical,
    covariance,
    domains,
    green,
    jets,
    recursion,
    riccati,
    transitions,
)

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"global-two-cone-Hadamard-state.json"
PARENT_SHA="4952891d95aa4c9f847b16b9771016c4c65a9eeb14cb932e9290ff03de530f9d"
sha,serialize=certificate.sha,certificate.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_proca_global_hadamard/*.py"))+sorted(ROOT.glob("tests/*.py"))
        +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen global causal-support report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_98_fully_rebuilt":PARENT_SHA,
        "same_named_action_and_same_original_global_clock":True,
        "new_high_frequency_state_not_a_promotion_of_old_mu_covariance":True}


@cache
def build_report():
    prior=prior_checks()
    identities,gates=audit.residuals(),audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A global Hadamard-state domain, covariance or recursion gate failed")
    exact=certificate.certify_residuals(identities)
    return {"schema":1,"claim":"P8-S6.99.GLOBAL_TWO_CONE_HADAMARD_STATE","date":"2026-09-09",
        "status":"COMPLETE_GLOBAL_REDUCED_QUADRATIC_SCALAR_THEORY_ON_PRESERVED_RETUNED_BOUNCE_ADMITS_A_POSITIVE_TWO_CONE_HADAMARD_STATE_WITH_EXACT_CAUSAL_CCR; RENORMALIZED_PHYSICAL_STRESS_INTERACTIONS_V_G_B_AND_ORIGINAL_P8_OPEN",
        "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
        "formulation":"FORMULATION.md",
        "written_proofs":["notes/green.md","notes/canonical.md","notes/riccati.md","notes/state.md","notes/wavefront.md","notes/scope.md"],
        "exact_residuals":exact,"named_exact_check_count":len(identities),
        "checked_scalar_entries":sum(value.rows*value.cols if isinstance(value,sp.MatrixBase) else 1
                                     for value in identities.values()),
        "proof_checks":gates,
        "physical_volume_operator_and_Green_CCR_conventions":serialize(green.data()),
        "fresh_complete_global_canonical_Laurent_chart":serialize(canonical.data()),
        "actual_initial_slab_graph_and_center_jets":serialize(jets.data()),
        "all_order_Riccati_inverse_and_four_order_native_tests":serialize({
            "principal":riccati.principal(),"coupled_native_recursion":recursion.data(),
            "global_mode_recursion":recursion.mode_data()}),
        "exact_positive_graph_Cauchy_covariance_and_infrared_completion":serialize(covariance.data()),
        "explicit_all_finite_strip_frequency_gaps":serialize(domains.data()),
        "complete_order_zero_chart_transitions":serialize(transitions.data()),
        "global_two_cone_wavefront_condition":{
            "V_plus":"xi!=0 and tau=-c(u)*|xi|/a(u) or tau=-|xi|/a(u)",
            "V_minus":"-V_plus",
            "positive_type_two_point_WF":"WF(W) subset V_plus x V_minus",
            "CCR_implies_decomposability":"WF(E_Q) subset (V_plus x V_minus) union (V_minus x V_plus)",
            "actual_global_operator":"Q=a^-3 Omega (partial_u-M_density)",
            "short_distance_proof":"all-orders Riccati symbol, Borel sum, exact Duhamel smoothing, signed diagonal phases and complete order-zero overlap projector uniqueness",
            "not_an_automatic_normally_hyperbolic_existence_theorem":True,
            "not_a_single_physical_metric_Klein_Gordon_Hadamard_condition":True},
        "controls":audit.controls(),
        "verdict":"The complete reduced scalar phase on the preserved global epsilon=1/200 bounce admits a fresh positive quasifree state with exact physical-volume CCR and the generalized two-cone Hadamard wavefront condition. The actual density operator has causal advanced/retarded Green operators. A new canonical Laurent bridge gives a positive negative-frequency graph; its all-order symmetric Riccati recursion, Borel sum and convex infrared completion define a positive state evolved by the exact global equation. Finite-strip frequency gaps and full order-zero chart overlaps preserve the signed microlocal subspaces. This is a written all-order state theorem with native exact algebraic and domain inputs, not a finite-order adiabatic approximation, ordinary free-matter stress calculation, interacting or UV theorem.",
        "not_established":["A preferred unique vacuum, numerical Borel threshold, global Fock unitary history or ground state at every momentum",
            "Ordinary one-metric Klein-Gordon Hadamard state, covariant renormalized physical stress, Ward identities or self-consistent semiclassical backreaction",
            "Interacting microcausality, canonical interaction/loop/omitted-operator/threshold control, or a quantitative quantum EFT validity band",
            "Full quantum tensor/vector system, physical vacuum amplitudes, finite-gravity dispersion or common-parent V/G/B matching",
            "Global continuation of the off-clock nearby solution, original realistic-field P8(a), original UV P8(b) classification or original P8 closure"],
        "verification_boundary":"Native literal global density operator and physical-volume Green sign, exact scaled symplectic and full Laurent identities, actual leading kinetic/gradient bridge and nonzero center Riccati derivative jets, exact universal frequency-sum Sylvester inversion and four-order coupled recurrence, noncommuting positive Gaussian Gram/CCR identities, explicit all-finite-strip positive-polynomial gap bounds, and complete rational order-zero chart maps with inverses. Pinned written Green, all-integer Riccati induction, parameter-dependent Borel construction, Duhamel smoothing, infrared state and global signed-wavefront/projector proofs. Every ancestor, report field and source/proof/test byte is independently replayed; scientific SymPy is unchanged. Not a machine-formalized proof of microlocal analysis."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The global two-cone Hadamard-state report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.99.GLOBAL_TWO_CONE_HADAMARD_STATE replay passed; physical stress, interactions, V/G/B and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
