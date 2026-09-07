"""Read-only source-preserving zero-relative tree-branch matching replay."""

import argparse
import hashlib
import json
from fractions import Fraction as Q
from functools import cache
from pathlib import Path

import sympy as sp
from p8 import verify_all as classification
from p8_composite_response import verify as prior

from . import branch, cayley, independent, matching, vacuum

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"composite-vacuum.json"
PRIOR_SHA = "34fa2ee35c2b33557e84b8fc181a202bed4793ed00d6aa8b76964ebbdc5dbe12"
CLASSIFICATION_SHA = "4e82b45d3d1daed0a5e4698f1dd39fa51e2eb09b6efcf3cf5caac3a12d45dec8"
CD_M1_SHA = "caf8c8e688a7565b9d00f921c099a28da00f97522ed26ad182a8227eb80cd4dd"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@cache
def prior_checks():
    reports = classification.ROOT/"certificates"
    pins = ((prior.REPORT, PRIOR_SHA),
            (reports/"classification.json", CLASSIFICATION_SHA),
            (reports/"witness-CD_matter.json", CD_M1_SHA))
    for path, expected in pins:
        if sha(path) != expected:
            raise ValueError(f"Pinned input changed: {path.name}")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    classification.check_reports()
    return {"S6_11_composite_response": PRIOR_SHA,
            "original_32_row_classification": CLASSIFICATION_SHA,
            "original_CD_matter_witness": CD_M1_SHA}


def controls():
    c, v, b, m = cayley.controls(), vacuum.controls(), branch.controls(), matching.controls()
    M2, m2 = cayley.M2, cayley.m2
    eps = next(symbol for symbol in v["pinned_bounce_initial_potential_over_M2m2"].free_symbols
               if symbol.name == "epsilon")
    initial_v = v["pinned_bounce_initial_potential_over_M2m2"]
    # V_b >= -epsilon, since each term on the right is nonnegative.
    positive_remainder = sp.Rational(9, 2)*eps**2+eps**2*(2+eps)/(1+eps)**2
    if sp.cancel(initial_v+eps-positive_remainder) != 0:
        raise ValueError("Same-potential lower-bound identity failed")
    polynomial = sp.Poly(v["constant_vacuum_potential_mismatch_cleared"], eps)
    if polynomial.all_coeffs() != [36, 72, 39, -2, 3]:
        raise ValueError("The constant-vacuum potential mismatch disappeared")
    domain_rejected = False
    try:
        matching.null_residual_budget(0, Q(1, 5))
    except ValueError:
        domain_rejected = True
    values = {
        "asymmetric_beta0_relative_tadpole_over_M2m2": sp.cancel(c["asymmetric_beta0_linear_relative_tadpole"]/(M2*m2)),
        "singular_Cayley_plus_determinant": c["singular_Cayley_plus_example"],
        "negative_root_outside_positive_chart": c["negative_root_not_positive_chart"],
        "bare_beta2_wrong_physical_mass_over_m2": sp.cancel(v["bare_beta2_wrong_physical_mass"]/m2),
        "omitted_vacuum_shift_mass_error_over_m2": sp.cancel(v["lost_vacuum_shift_mass_error"]/m2),
        "wrong_clock_mass_over_m2": sp.cancel(v["wrong_equal_time_clock_mass"]/m2),
        "positive_mass_Dirichlet_resonant_equation": b["positive_mass_Dirichlet_resonant_equation"],
        "resonant_left_boundary": b["positive_mass_Dirichlet_left"],
        "resonant_right_boundary": b["positive_mass_Dirichlet_right"],
        "nonzero_solution_with_zero_Dirichlet_data": b["positive_mass_Dirichlet_nonzero_solution"],
        "nonzero_heavy_initial_state": b["nonzero_heavy_initial_state"],
        "branch_switch_inverse_loss": b["branch_switch_loses_linear_inverse"],
        "symmetry_breaking_relative_source": b["broken_symmetry_linear_source"],
        "even_heavy_operator_can_generate_light_loop_dependence": b["even_heavy_quadratic_coefficient_can_produce_light_loop_dependence"],
        "conditional_not_computed_CD_residual_fixture": b["conditional_residual_fixture"],
        "conditional_ball_fixture": b["conditional_ball_fixture"],
        "pure_heavy_degree_2_component_minimum_loops": branch.closed_heavy_graph_loop_bound([2, 2, 2]),
        "heavy_degree_2_and_4_component_minimum_loops": branch.closed_heavy_graph_loop_bound([2, 4]),
        "same_potential_gap_for_epsilon_at_most_1_over_10000": m["constant_vacuum_vs_pinned_bounce_potential_minimum_gap"],
        "dropping_exact_free_chi_understates_null_budget": m["dropping_free_chi_understates_null_residual_budget"],
        "static_GR_sharp_endpoint_error_threshold": m["static_GR_endpoint_error_threshold"],
        "NEC_violation_is_outside_null_obstruction": m["NEC_violating_source_can_reverse_GR_Hdot"],
        "exact_actual_null_residual_budget": matching.null_residual_budget(0, 0),
        "approximate_actual_null_residual_fixture": matching.null_residual_budget(Q(1, 10), Q(1, 100)),
        "oversized_chi_error_rejected": domain_rejected,
    }
    expected = (sp.Rational(-1, 16), 0, -1, 1, sp.Rational(3, 4), 1,
                0, 0, 0, 1, 1, 0, 1, 1, Q(2, 25), True, 1, 2,
                Q(3749, 10000), Q(1, 100), Q(8, 5), -1, Q(801, 100), Q(78081, 10000), True)
    if tuple(values.values()) != expected:
        raise ValueError("A parity, source, vacuum, functional-inverse or matching control failed")
    return {"checked_values": {key: str(value) for key, value in values.items()},
            "same_potential_gap_positive_remainder_identity": "0",
            "same_potential_mismatch_polynomial_coefficients_descending": [36, 72, 39, -2, 3]}


@cache
def build_report():
    previous = prior_checks()
    residuals = {"exact_Cayley_and_full_potential": cayley.checks(),
                 "physical_vacuum_mass_source_and_shifts": vacuum.checks(),
                 "conserved_physical_source_projector": vacuum.source_projector_checks(),
                 "stationary_branch_and_complete_differential_operator": branch.checks(),
                 "original_CD_M1_physical_and_operator_dictionary": matching.checks()}
    if any(sp.simplify(value) != 0 for group in residuals.values() for value in group.values()):
        raise ValueError("An exact source, branch or matching identity failed")
    rational = independent.checks()
    omissions = controls()
    v = vacuum.derive()
    positives = {"physical_Planck_squared_over_M2": sp.cancel(v["physical_Planck_squared"]/v["M2"]),
                 "relative_Planck_squared_over_M2": sp.cancel(v["relative_Planck_squared"]/v["M2"]),
                 "physical_FP_mass_squared_over_m2": sp.cancel(v["physical_FP_mass_squared"]/v["m2"]),
                 "relative_auxiliary_probe_residue_times_M2": sp.cancel(v["relative_auxiliary_probe_residue"]*v["M2"])}
    if any(value <= 0 for value in positives.values()):
        raise ValueError("The conditional flat-vacuum mass or kinetic sign failed")
    sources = sorted(ROOT.glob("src/p8_composite_vacuum/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1,
        "claim": "P8-S6.12.COMPOSITE",
        "date": "2026-09-06",
        "status": "EXACT_SOURCE_PRESERVING_ZERO_RELATIVE_CANONICAL_TREE_BRANCH_MATCHING_OBSTRUCTION; OTHER_PARENT_BRANCHES_AND_P8_OPEN",
        "prior_report_sha256": previous,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
        "formulation": "FORMULATION.md",
        "written_proof": "notes/proof.md",
        "source_audit": "notes/sources.md",
        "prior_art": {"primary_source": "https://arxiv.org/pdf/1409.3146",
                      "audited_equations": "2.7,4.11-4.15,4.20-4.22,4.23-4.35",
                      "boundary": "Proportional GR/source results are prior art; exact selected-branch proof is not generic nonlinear equivalence or a novelty claim"},
        "exact_residuals": {name: dict.fromkeys(group, "0") for name, group in residuals.items()},
        "independent_Fraction_replay": rational,
        "negative_and_positive_controls": omissions,
        "conditional_flat_vacuum_positive_coefficients": {key: str(value) for key, value in positives.items()},
        "model_and_frame": {
            "gravity": "Equal positive constants M², interaction scale m>0, bare beta=(0,0,1,0,0), +--- convention; each EH term -M² sqrt|g_i| R[g_i]/2",
            "physical_metric": "mathcalG=g_eff=g(I+sqrt(g^-1 f))² with fixed alpha=beta=1; the literal full light and external source actions depend on mathcalG alone",
            "chart": "Delta=(I-S)(I+S)^-1 is mathcalG-self-adjoint; I+-Delta invertible, g/f Lorentzian and selected real positive-root chart connected to Delta=0",
            "map": "g=mathcalG(I+Delta)²/4, f=mathcalG(I-Delta)²/4, S=(I+Delta)^-1(I-Delta)",
            "exchange": "g<->f sends Delta->-Delta and leaves mathcalG, every prescribed light field and physical source fixed",
            "general_parity_scope": "Same stationary identity for constant beta_n=beta_(4-n), equal Einstein coefficients, exchange-compatible gravitational boundary terms and physical-metric-only sources",
        },
        "distinct_light_classes": {
            "pinned_parent": "One shared positive-kinetic canonical scalar with its reconstructed rolling potential; no automatic identification with old CD/M1",
            "explicit_M1_extension": "Separately named class with positive field metric Z_AB(psi), V(psi), plus exact independent free chi kinetic term (partial chi)²/2; no chi potential, mixed kinetic or chi dependence in Z,V",
            "tree_input_restriction": "Canonical light action contains no preinserted C/D or higher-derivative operators; arbitrary prescribed inputs would instead be retained",
            "potential_matching": "Any vacuum and rolling claims must concern one specified action/potential; no off-tube extension or same-potential construction is supplied here",
        },
        "exact_tree_branch": {
            "potential": "-M²m² sqrt|mathcalG| [6-2e2(Delta)+6e4(Delta)]/16",
            "relative_equation": "E_Delta[mathcalG,0,psi;J]=0 off shell in all light histories and sources, retaining metric/time derivatives",
            "action": "S[mathcalG,0,psi;J]=integral sqrt|mathcalG|[-M²R[mathcalG]/4-3M²m²/8+L_light]+Sprobe[mathcalG;J]",
            "stationarity": "Light variation of the substituted action equals the full light variation because the functional chain-rule term multiplies E_Delta=0",
            "no_new_tree_operators": "No relative-sector-induced C/D DHOST or curvature-squared operators at any derivative order on this exact branch; existing light GR/matter trees remain",
            "diagrammatic_check": "Every vertex has even relative degree; a nonempty connected heavy subgraph with no heavy external legs has E>=V and at least one loop",
            "loops": "Evenness permits light-dependent heavy log determinants; no quantum cancellation or radiative closure is established",
        },
        "functional_elimination_contract": {
            "causal_option": "Specify regular gauge/constraints and a well-posed relative initial-value problem with zero relative canonical and auxiliary data; uniqueness then selects Delta=0",
            "stationary_option": "A continuously differentiable relative equation on specified Banach spaces, continuous selected solution functional, common zero-relative boundary/state prescription, and an isomorphism of its full constrained linearization along a connected parameter domain",
            "continuation": "The selected zero branch is open by local uniqueness and closed by continuity, hence persists on the connected domain; this is configuration-space continuation, not a required physical-time vacuum-to-bounce trajectory",
            "not_supplied": "An actual rolling-CD constrained inverse norm, nonlinear well-posedness, state selection, or same-potential off-tube extension",
            "conditional_residual": "If E(h)=Lh+N(h)=r, N(0)=0, ||L^-1||<=K and ||L^-1(N(h)-N(k))||<=eta||h-k|| with eta<1 in a complete ball, then ||h||<=K||r||/(1-eta)",
            "conditional_existence": "K||r||<=(1-eta)*radius makes that contraction map the ball to itself; residual norm includes any boundary/initial-data mismatch",
            "fixture_scope": "The exact API's rational K=2,eta=1/4,r=3/100 calculation is illustrative arithmetic, not a computed bimetric-CD estimate",
            "failure_controls": "Positive algebraic mass does not exclude resonant boundary kernels, nonzero heavy data, distinct saddles, or an inverse loss at branch joining",
        },
        "conditional_proportional_vacuum": {
            "required_fixed_potential_jet": "V(psi_v)=-3M²m²/8, V_A(psi_v)=0 and canonically normalized Hessian nonnegative; positive clock mass when that vacuum contract requires it, free chi remains massless",
            "all_beta_vacuum_shift": "beta_eff=(-3/8,-3/8,5/8,-3/8,-3/8)",
            "physical_geometry_and_clock": "g=f=mathcalG/4, physical mathcalG Minkowski, T=2T_g; conventional g-clock relative mass m² becomes physical m²/4",
            "full_FP_invariant": "-M²m²[tr(Delta²)-(tr Delta)²]/16, with five constrained source-free massive polarizations at flat quadratic order",
            "physical_source_action": "For eij eij=1, delta mathcalG_ij=-gamma eij and delta Sprobe=(1/2) integral T^ij delta mathcalG_ij, the term j gamma/4 has j=-2Pi_TT",
            "TT_quadratic_action": "-[M_eff² D gamma²+M_relative²(D+m_FP²)Delta_TT²]/8+j gamma/4; M_eff²=M²/2,M_relative²=2M²,m_FP²=m²/4,D=partial_T²+kphys²",
            "physical_source_response": "gamma=j/(M_eff²D), relative physical-source residue zero; stationary source functional j²/(8M_eff²D)",
            "full_conserved_source": "Massless response carries (P2-P0/2)/(M_eff²D); massive P2 has zero projection onto physical sources here, not negative norm or zero mass",
            "independent_relative_probe": "A diagnostic nonphysical relative source would see positive residue 1/(2M²) and the massive pole; it is not added as a matter coupling",
            "same_potential_negative_control": "Pinned rolling centre V_b/(M²m²)=9epsilon²/2-epsilon/(1+epsilon)² differs from constant vacuum -3/8 by at least 3749/10000 for 0<epsilon<=1/10000",
            "scope": "Flat quadratic vacuum audit under these potential-jet assumptions; not existence of the required same-potential vacuum and rolling action, nor full-parent health",
        },
        "original_CD_M1_matching_obstruction": {
            "operator_dictionary": "Frozen original physical metric and nondifferential clock dictionary, original M=tau=1 coefficient normalization; target F2_X(0,1)=-1/2,A3(0,1)=1 whereas this canonical zero branch has both zero",
            "operator_remainder_budgets": "Target=branch+remainder requires |R_F2X(0,1)|>=1/2 and |R_A3(0,1)|>=1; not invariant under unaccounted derivative or matter-metric changes",
            "physical_normalization": "M_target²=M_eff²=M²/2, H_CD=4T/(tau²+T²), chidot_CD=M_target/[10tau(1+(T/tau)²)^6]",
            "NEC_null_equation": "-2M_target² Hdot=Z_AB psidot^A psidot^B+chidot²>=0 on this canonical branch",
            "finite_window_metric_error": "On [-tau/2,tau/2], any Hubble approximation with endpoint errors both strictly below 8/(5tau) would require increasing H and is impossible; no Hdot error premise needed",
            "sharp_endpoint_control": "Static GR H=0 reaches the endpoint threshold, but is not claimed to have the old nonzero rolling M1 charge",
            "actual_null_residual": "In -2M_target²Hdot=n_clock+chidot²+R_null with n_clock>=0, the exact target requires R_null<=-(801/100)M_target²/tau²",
            "approximate_residual_domain": "If |tau²Hdot-4|<=epsilon_Hdot<4 and |tau chidot/M_target-1/10|<=epsilon_chidot<1/10, both errors nonnegative, then -R_null/(M_target²/tau²)>=2(4-epsilon_Hdot)+(1/10-epsilon_chidot)²",
            "error_scope": "An actual omitted null-equation residual evaluated on the physical solution/field map; not a per-operator or loop smallness assertion without the needed norm dictionary",
        },
        "verification_boundary": "Source-hashed immutable replay, exact SymPy and independent coefficientwise Fraction algebra, full root covariant/source/FP audits, and explicit written functional theorems; not proof-assistant formalized",
        "not_established": [
            "An exclusion of nonzero-relative backgrounds, other boundary/data/state prescriptions or singular/noninvertible stationary branches",
            "An exclusion of the whole composite parent family, other parents, the old CD/M1 row or all ultraviolet possibilities",
            "A required physical-time history from vacuum to the bounce",
            "Actual CD heavy-operator invertibility, quantitative heavy-gap separation, cutoff hierarchy or nonlinear full-parent health",
            "Quantum or mixed-loop matching, positivity, radiative closure, or cancellation of allowed heavy loop operators",
            "A common reconstructed rolling potential and conventional vacuum, or silently adding the old free chi to an unchanged one-scalar background",
            "Completion of the adopted S6 matching gate, P8(b), P8(a), or P8",
        ],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Composite source-preserving zero-relative-branch certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.12.COMPOSITE: exact zero-relative canonical tree-branch matching replay passed; other branches and P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
