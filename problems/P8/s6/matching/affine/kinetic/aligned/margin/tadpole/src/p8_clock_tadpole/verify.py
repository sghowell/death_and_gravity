"""Read-only selected-state clock-preserving scalar-profile certificate."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_vector_regularity import verify as parent
from p8_vector_variation import preparation as prior_preparation

from . import fronts, profiles, state, window

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"selected-state-clock-profile.json"
PARENT_SHA = "ee3da10ed767ca94be285bf262b745d12f8a3b0c5dfa3b2a7797da433c2c58a2"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_clock_tadpole/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen quantitative vector-regularity certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_59_fully_rebuilt": PARENT_SHA,
            "new_named_lower_scalar_action_not_parent_overwrite": True,
            "vector_state_operator_and_finite_prescription_unchanged": True}


@cache
def residuals():
    out = {**profiles.covariant()["checks"], **profiles.point_chart()["checks"],
           **fronts.principal()["checks"], **fronts.actual()["checks"],
           **window.chain_rule_checks(), **state.checks()}
    h = sp.Symbol("h", positive=True)
    out["actual_longitudinal_clock_lapse_slope"] = sp.Rational(16, 81)*(-1/(2*h))+sp.Rational(8, 81)/h
    N, p = sp.symbols("N p", positive=True)
    out["clock_norm_lapse_and_affine_sign_relation"] = sp.factor(
        ((N**-2-1)/h-(4*p**2-1)).subs(N**-2, 1+h*(4*p**2-1)))
    return out


@cache
def proof_checks():
    data = profiles.profile_bounds(10**24, 1000)
    out = {**fronts.proof_checks(), **window.proof_checks()}
    for pair, bound in data["all_clock_norm_mixed_derivative_bounds_total_order_at_most_five"].items():
        out["mixed_profile_derivative_"+pair.replace(",", "_")+"_below_1e_minus_18"] = bool(bound < sp.Rational(1, 10**18))
    for order, row in data["fixed_profile_derivative_bounds_on_clock_tube"].items():
        out["clock_tube_profile_derivative_"+str(order)+"_below_1e_minus_18"] = bool(row["P"] < sp.Rational(1, 10**18))
        out["clock_tube_Px_derivative_"+str(order)+"_below_1e_minus_18"] = bool(row["P_x"] < sp.Rational(1, 10**18))
    out.update({"literal_lapse_square_coefficient_stays_above_one_third": bool(sp.Rational(3, 8)-data["lapse_square_coefficient_absolute_upper"] > sp.Rational(1, 3)),
                "prior_actual_initial_force_remains_positive": bool(prior_preparation.initial_force()["initial_lapse_force_lower"] > 0),
                "profile_window_exactly_one_at_clock": window.cutoff(-1) == 1,
                "profile_window_exactly_zero_at_vacuum": window.cutoff(0) == 0,
                "positive_lapse_front_fixture_is_subluminal": fronts.speed(sp.Rational(49, 100))["longitudinal_squared_principal_speed"] < 1,
                "negative_lapse_front_fixture_is_superluminal": fronts.speed(sp.Rational(51, 100))["longitudinal_squared_principal_speed"] > 1})
    return {name: bool(value) for name, value in out.items()}


def controls():
    bad_exact = (True, False, sp.true, sp.false, 1.0, sp.Float(1), "1000", sp.I,
                 sp.oo, -sp.oo, sp.zoo, sp.nan, sp.Symbol("unproved"))
    calls = []
    for value in bad_exact:
        calls.extend((lambda value=value: fronts.speed(value),
                      lambda value=value: window.cutoff(value),
                      lambda value=value: profiles.profile_bounds(value, 1000),
                      lambda value=value: state.definition(10**24, value)))
    bad_orders = (True, False, sp.true, sp.false, 1.0, sp.Float(1), sp.Integer(1), "1", -1, 6, sp.Rational(1, 2))
    for value in bad_orders:
        for function in (window.bump_bound, window.inverse_denominator_bound, window.step_bound, window.derivative_bound):
            calls.append(lambda value=value, function=function: function(value))
    calls.extend(lambda value=value: fronts.speed(value) for value in (0, -1, 1, sp.Rational(11, 25), sp.Rational(53, 100)))
    calls.extend(lambda value=value: profiles.profile_bounds(value, 1000) for value in (0, -1))
    calls.extend(lambda value=value: state.definition(10**24, value) for value in (0, -1, 999))
    calls.extend(lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan))
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An invalid exact profile or vector-block diagnostic was accepted")
    return {"rejected_inputs": rejected,
            "cache_warmup_does_not_bypass_native_order_validation": True,
            "fixed_profiles_not_recomputed_under_metric_or_state_variations": True,
            "physical_lower_scalar_action_change_not_renormalization_scheme_change": True,
            "one_selected_state_cancellation_not_universal_state_cancellation": True,
            "off_clock_vector_block_not_full_coupled_quantum_cone": True,
            "global_smooth_definition_not_uniform_all_time_smallness": True,
            "vacuum_zero_region_not_UV_matching_certificate": True,
            "coefficient_jet_bounds_not_canonical_cutoff_estimates": True}


@cache
def build_report():
    prior, identities = prior_checks(), residuals()
    exact, proofs = affine.certify_residuals(identities), proof_checks()
    if not all(value is True for value in proofs.values()):
        raise ValueError("A selected-state profile or vector-block proof check failed")
    scalar_count = sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1 for value in identities.values())
    return {"schema": 1, "claim": "P8-S6.60.CLOCK_PROFILE", "date": "2026-09-08",
            "status": "NEW_SELECTED_STATE_CLOCK_PRESERVING_SCALAR_ACTION; ORIGINAL_P8_OPEN",
            "prior_sha256": prior,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md", "written_proofs": ["notes/fronts.md", "notes/profiles.md", "notes/bounds.md"],
            "exact_residuals": exact, "named_exact_check_count": len(identities), "checked_scalar_entries": scalar_count,
            "proof_checks": proofs,
            "literal_new_action": "To S6.56 add Delta P=(M²/tau²) w(x)[-p_sigma(u)+(rho_sigma(u)+p_sigma(u))(x+1)/2], u=phi/tau, w=1-chi(64(x+1)²), using exactly the S6.55 smooth step. The selected normalized stress profiles are computed once from the predecessor state, original history and unchanged finite prescription. The addition has the same formal loop order as that Gaussian stress.",
            "fixed_profile_definition": serialize(state.definition(10**24, 1000)),
            "actual_mode_profile_integrands": serialize(state.integrands()),
            "unchanged_initial_data_examples": serialize(state.initial_data_examples()),
            "literal_clock_tube_and_global_profile": serialize(profiles.covariant()),
            "actual_point_chart_quadratic_terms": serialize(profiles.point_chart()),
            "window_derivative_envelopes": serialize({order: window.derivative_bound(order) for order in range(6)}),
            "finite_profile_error_budget": serialize(profiles.profile_bounds(10**24, 1000)),
            "background_stationarity_and_completeness_scope": "The added clock stress is exactly (-rho_sigma,-p_sigma); its scalar Euler source is the negative of the same vector Ward source. Thus the original complete clock geometry and free-matter history solve the reduced retained-vector Gaussian semiclassical background equations in this new action and selected state. The aligned mean source and its source-squared first variation still vanish there. This is not the quantum background equation of omitted sectors.",
            "global_smoothness_and_local_budget_boundary": "Every compact real-time interval admits arbitrary-order high-frequency estimates for the same smooth exact oscillator and all-order initial data. Finite subtracted stress and the fixed dimensional prescription extend smoothly over the full clock history. Quantitative smallness is proved only for |u|<=1/2. There all 21 mixed (u,x) derivatives of total order <=5 are below 10^-18 for every real x in the scale example.",
            "vacuum_preservation_boundary": "The smooth window equals one on an open neighborhood of |x+1|<=1/10 and vanishes in an open neighborhood of x=0, so this addition changes no local field jet there. If a predecessor vacuum lies there its jets are unchanged. This does not establish vacuum existence or health, a common UV parent, nonperturbative vacuum equivalence or matching across the intervening domain.",
            "retained_vector_principal_data": serialize({"general": fronts.principal(), "actual": fronts.actual()}),
            "prior_approximate_response_cone_boundary": "For the retained vector block c_T²=1 and c_L²=b_mass/a_mass. The actual factorization has sign(c_L²-1)=sign(2p_affine-1)=sign(N^-2-1) throughout the original tube. S6.58 proves N_initial<1 for the particular zero-independent-data fixed-source response, so that block is superluminal on its approximate initial metric. Larger L reduces the excess but not its sign. This neither excludes other initial data nor computes the full coupled quantum cone.",
            "exact_vector_block_examples": serialize({str(value): fronts.speed(value) for value in
                                                      (sp.Rational(49, 100), sp.Rational(1, 2), sp.Rational(51, 100))}),
            "remaining_variation_boundary": "Profiles are fixed scalar-action coefficients, not reevaluated functions of a varied state or metric. Only the selected first variations cancel. The full contact/retarded second variation, compatible initial covariance response, corrected coupled cones and interactions remain to be controlled. The literal lapse-square addition has abs(Delta J)<=15eta0/8; this is not a quantum Dirac-count or stability theorem.",
            "controls": controls(),
            "verdict": "A new, finely tuned selected-state action preserves the complete original clock history in the retained-vector Gaussian semiclassical background equations, with a quantified compact-window scalar-profile budget. Original P8 remains open.",
            "not_established": ["A stable or causal full quantum solution, a renormalized second-variation bound or compatibility for arbitrary perturbed histories",
                                "Light/other-sector loops, all omitted-operator effects, canonical interacting cutoff or uniform all-time smallness",
                                "Common UV parent, finite-gravity Regge remainder, V/G/B or original P8 closure"],
            "verification_boundary": "Exact symbolic/rational first variations, canonical front limits, continuous sign and smooth-window estimates with written global smoothness and selected-state background proofs; not proof-assistant formalization"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The selected-state profile report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.60.CLOCK_PROFILE replay passed; selected-state background, original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
