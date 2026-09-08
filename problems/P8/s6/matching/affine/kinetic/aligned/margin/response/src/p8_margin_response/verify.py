"""Read-only fixed-vector-stress homogeneous response certificate."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_aligned_margin import verify as parent

from . import bounds, system

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"fixed-vector-stress-response.json"
PARENT_SHA = "b6952b45787005033a7ae19176c4668525044b37a2afb8753651ca99915d47dd"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_margin_response/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen clock-cone margin certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_56_fully_rebuilt": PARENT_SHA,
            "action_state_and_finite_matching_unchanged": True,
            "frozen_vector_first_variation_bounds_not_quantum_feedback_norm": True}


def residuals():
    return {**system.checks(), **bounds.checks()}


def controls():
    bad = (True, False, sp.true, sp.false, 1.0, sp.Float(1), "1", sp.I,
           sp.oo, -sp.oo, sp.zoo, sp.nan, sp.Symbol("unproved"), -1, sp.Rational(1, 10**5))
    calls = [lambda value=value: bounds.response_bounds(value, 0) for value in bad]
    calls += [lambda value=value: bounds.response_bounds(0, value) for value in bad]
    calls += [lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan)]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An invalid exact source budget or report value was accepted")
    return {"rejected_inputs": rejected,
            "homogeneous_action_restricted_before_nonzero_q_shift_elimination": True,
            "physical_metric_lapse_variation_included": True,
            "auxiliary_lapse_not_an_independent_zero_initial_datum": True,
            "joint_free_matter_momentum_preserved": True,
            "positive_chart_lift_not_a_nonlinear_solution": True,
            "fixed_source_not_a_self_consistent_quantum_feedback_map": True,
            "no_new_counterterm_or_discarded_causal_branch": True}


@cache
def build_report():
    pins, identities = prior_checks(), residuals()
    exact, proofs = affine.certify_residuals(identities), bounds.proof_checks()
    if not all(value is True for value in proofs.values()):
        raise ValueError("A continuous response or physical geometry bound failed")
    scalar_count = sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1 for value in identities.values())
    return {"schema": 1, "claim": "P8-S6.57.FIXED_STRESS_RESPONSE", "date": "2026-09-08",
            "status": "CONTROLLED_FIXED_VECTOR_STRESS_FIRST_ORDER_HOMOGENEOUS_RESPONSE; ORIGINAL_P8_OPEN",
            "prior_sha256": pins,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md",
            "written_proofs": ["notes/response.md", "notes/geometry.md", "notes/literature.md"],
            "exact_residuals": exact, "named_exact_check_count": len(identities),
            "checked_scalar_entries": scalar_count, "proof_checks": proofs,
            "literal_model_state_and_source": "Keep S6.56 with 0<epsilon_margin<=1/100 and the unchanged S6.55 Hadamard vector state and finite prescription. On I=[-1/2,1/2], use its actual fixed-background vector energy and pressure as first-order homogeneous forces. Normalize values by M²/tau² and first physical time derivatives by M²/tau³.",
            "physical_source_variation": "The literal homogeneous action is restricted before any nonzero spatial-mode shift division. Physical curvature is v+delta*n, delta=1/(2h). The stress variation is F_v*v+F_n*n with F_v=3p,F_n=3delta*p-rho. The compatible clock source is the inherited J_clock=-[rho'+3H(rho+p)], not an independent adjustment.",
            "regular_canonical_reduction": "Direct Legendre transformation gives H=-(p_v+3ell*s)²/12+p_m²/2+[Theta(p_v+3ell*s)-w*p_m-F_n]²/[4(J+delta_J)]-F_v*v. There is no Theta or spatial inverse. Zero independent canonical initial perturbations leave n(u0)=-F_n(u0)/[2(J+delta_J)]. The exact invariant a³(p_m+3ell*v) remains zero.",
            "normalized_two_component_system": "For S=p_v+3ell*s,Y=4a³S, the system for (v,Y) has matrix [[B,A/(4a³)],[-4a³C,-B]] and forcing [-Theta*F_n/(2J_e),4a³(F_v+3ell*w*F_n/(2J_e))], with A=Theta²/(2J_e)-1/6,B=3Theta*w*ell/(2J_e),C=9ell²[1+w²/(2J_e)]. All entries are derived from the changed literal action.",
            "continuous_inverse_bound": "J_e>3/8,1<=a³<4,abs(Theta)<=2,abs(ell)<=1/10,abs(w)<=1/20. The two matrix row bounds are 283/200 and 928/625, both below 3/2. For value budget eta0, forcing is below 49eta0. Exact exp-series control gives exp(3/2)<5, hence the zero-data phase response is <=131eta0 and its derivative <=246eta0 over the full interval.",
            "reconstruction_and_time_derivative_bounds": "The actual algebraic lapse and preserved matter momentum give abs(n)<=94eta0,abs(s)<=44eta0. Exact background derivative envelopes give abs(J_e')<185 and a lapse-coefficient derivative sum below 346. With first source derivative budget eta1, abs(n')<=48000eta0+4eta1 and abs(v')<=193eta0.",
            "exact_positive_physical_chart_lift": "Set N_app=1+n,a_hat,app=a exp(v),a_phys,app=a exp[v+omega(u,N_app)] with the original omega=-log[(h-1+N^-2)/h]/4. For 0<=eta0,eta1<=10^-6 the lapse stays inside the original X tube. The physical fractional scale change is <=450eta0 and tau times the actual Hubble difference is <=98000eta0+8eta1.",
            "approximate_metric_turning_minimum": "The lifted metric retains negative Hubble rate at u=-1/4 and positive Hubble rate at u=1/4 because the old magnitudes are 16/17 and the permitted correction is below 1/10. It therefore has an interior local scale-factor minimum. This concerns the specified approximate metric, not an exact quantum solution or a proved nondegenerate/unique bounce.",
            "exact_derivative_majorants": serialize(bounds.derivative_bounds()),
            "exact_operator_and_reconstruction_constants": serialize(bounds.constants()),
            "rounded_physical_scale_example": serialize(bounds.response_bounds(sp.Rational(1, 10**14), sp.Rational(1, 10**14))),
            "actual_vector_scale_example": serialize(bounds.vector_example(10**12, 1000)),
            "quantum_feedback_boundary": "Small fixed-background sources do not bound the causal second functional variation on perturbed histories. The local contact and retarded commutator contributions, compatible Hadamard preparation, constraints and nonlinear residual still require control. No quantum feedback contraction or full semiclassical residual is claimed.",
            "literature_scope": "Anderson, Molina-Paris and Mottola, arXiv:gr-qc/0209075, section III, distinguishes local contact and retarded response terms. Meda and Pinamonti, arXiv:2201.10288v2, treats a scalar toy system, not this Proca bounce. Neither is imported as a numerical response or stability theorem here; exact links and boundaries are in notes/literature.md.",
            "controls": controls(),
            "verdict": "The actual bounded vector source has a unique controlled first-order homogeneous response through the bounce, and its specified physical chart lift remains close and has a turning minimum. Full quantum feedback and original P8 remain open.",
            "not_established": ["A self-consistent semiclassical solution or quantitative nonlinear/quantum residual on the lifted metric",
                                "The complete causal second functional variation and corrected perturbation cones or interacting cutoff",
                                "A global quantum-complete history, unique/nondegenerate corrected bounce, V/G/B or original P8 closure"],
            "verification_boundary": "Exact symbolic/rational controls and written continuous linear-ODE, algebraic reconstruction and physical point-chart estimates; not proof-assistant formalization or a full quantum solution"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The fixed-vector-stress response report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.57.FIXED_STRESS_RESPONSE replay passed; fixed-source linear response, original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
