"""Read-only controlled Hadamard vector Cauchy-state certificate."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_vector_evolution import verify as parent

from . import cutoffs, preparation, proca, series, transfer

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"controlled-hadamard-vector-state.json"
PARENT_SHA = "845f1ca80239f75a03d59d5a4cfcb05b3e2fdeefe084e49b71a57480e7087a47"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_vector_hadamard/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen differentiated vector clock-source certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_54_fully_rebuilt": PARENT_SHA,
            "frozen_fourth_order_state_not_relabelled": True,
            "original_action_physical_clock_and_finite_counterterms_unchanged": True}


def residuals():
    out = {**series.low_order_checks(), **cutoffs.checks(), **cutoffs.structural_checks(),
           **proca.minkowski_seed(), **proca.reconstruction_checks()}
    for kind in ("transverse", "longitudinal"):
        out.update({kind+"_order_eight_"+key: value for key, value in series.coefficient_bounds(kind, 4).items()
                    if key.endswith("reconstruction")})
    return out


def proof_checks():
    out = {**cutoffs.proof_checks(), **transfer.proof_checks(), **transfer.scale_checks()}
    for kind in ("transverse", "longitudinal"):
        data = series.coefficient_bounds(kind, 4)
        threshold = cutoffs.threshold(kind, 4, 1000)
        out[kind+"_first_higher_threshold_doubles"] = bool(threshold >= 2*cutoffs.threshold(kind, 3, 1000))
        out[kind+"_first_higher_threshold_controls_both_slopes"] = bool(threshold >= 16*(1+max(
            data["coefficient_upper"], data["frequency_coefficient_slope_upper"])))
    return out


def controls():
    bad = (True, False, sp.true, sp.false, 1.0, sp.Float(1), "1", sp.I,
           sp.oo, -sp.oo, sp.zoo, sp.nan, sp.Symbol("unproved"), 0, -1)
    calls = [lambda value=value: transfer.physical_bounds(value, 1000) for value in bad]
    calls += [lambda value=value: transfer.physical_bounds(10**12, value) for value in (*bad, 999)]
    bad_orders = (True, False, sp.true, sp.false, 1.0, sp.Float(1), "1", -1, sp.Integer(1), sp.nan, sp.I, None)
    calls += [lambda value=value: series.coefficient("transverse", value) for value in bad_orders]
    calls += [lambda value=value: cutoffs.bump_derivative_polynomial(value) for value in bad_orders]
    calls += [lambda value=value: cutoffs.threshold("transverse", value, 1000) for value in
              (True, 1.0, sp.Integer(4), 0, 1, 2)]
    calls += [lambda value=value: cutoffs.threshold("transverse", 3, value) for value in (*bad, 999)]
    calls += [lambda value=value: cutoffs.turn_on(value) for value in bad[:13]]
    calls += [lambda value=value: preparation.initial_data("transverse", 1000, value) for value in (*bad, 999)]
    calls += [lambda value=value: series.coefficient(value, 0) for value in ("scalar", None, True)]
    calls += [lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan)]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An invalid state scale, order, cutoff or report value was accepted")
    return {"rejected_inputs": rejected,
            "locally_finite_definition_not_fixed_order_numerical_truncation": True,
            "initial_frequency_and_slope_both_controlled": True,
            "cutoffs_use_time_independent_comoving_frequency": True,
            "three_positive_constrained_Proca_polarizations_not_four_scalars": True,
            "auxiliary_metric_is_reference_proof_device_not_physical_bounce_replacement": True,
            "Hadamard_reference_constants_not_used_as_quantitative_error_coefficients": True,
            "finite_state_change_adds_no_pole_or_counterterm": True}


@cache
def build_report():
    pins, identities = prior_checks(), residuals()
    exact, proofs = affine.certify_residuals(identities), proof_checks()
    if not all(value is True for value in proofs.values()):
        raise ValueError("A cutoff, state-transfer or regulator bound failed")
    scalar_count = sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1 for value in identities.values())
    return {"schema": 1, "claim": "P8-S6.55.VECTOR_HADAMARD", "date": "2026-09-08",
            "status": "NEW_CONTROLLED_HADAMARD_VECTOR_STATE_WITH_MATCHED_FIRST_VARIATION_BOUNDS; ORIGINAL_P8_OPEN",
            "prior_sha256": pins,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md",
            "written_proofs": ["notes/construction.md", "notes/hadamard.md", "notes/transfer.md"],
            "exact_residuals": exact, "named_exact_check_count": len(identities),
            "checked_scalar_entries": scalar_count, "proof_checks": proofs,
            "literal_component": "A new cutoff-summed all-order Cauchy state of the unchanged retained free Proca vector on the actual clock. Physical frame, action and S6.53 finite matching prescription are unchanged. The old S6.50 fourth-order Gaussian is not relabelled as Hadamard.",
            "all_order_recurrence": "For S=1+sum P_n t^n, t=omega^-2, L=lambda+D S/S, cancel rho_res=(1-S²)/t-U-D L/2+L²/4. At order t^(n-1) the new coefficient enters only as -2P_n. Formal inverse and convolution recurrences define each P_n exactly. The actual rational coefficient ring is preserved under the fixed-comoving derivative, giving finite computable coefficient and slope envelopes at every finite order.",
            "locally_finite_state_definition": "nu²=m²+k²/(25/16)² is time independent. Lambda_3=2m; Lambda_n=max(2Lambda_(n-1),2^n[1+C_n]), where C_n is the larger exact box bound on P_n and D0 P_n+(1-2n)lambda P_n. Add chi(nu/Lambda_n) P_n omega^(1-2n) for every n>=3 to W_4. Chi is the explicit smooth flat turn-on from zero at one to one at two. Threshold doubling makes each momentum sum finite; the same rule proves all-order value and slope asymptotics.",
            "sixth_order_and_full_correction_envelopes": serialize({kind: series.leading_correction_bounds(kind) for kind in ("transverse", "longitudinal")}),
            "first_higher_cutoff_controls": serialize({kind: {"order_eight_bounds": series.coefficient_bounds(kind, 4),
                 "Lambda_3_at_mass_1000": cutoffs.threshold(kind, 3, 1000),
                 "Lambda_4_at_mass_1000": cutoffs.threshold(kind, 4, 1000)} for kind in ("transverse", "longitudinal")}),
            "explicit_Cauchy_data_examples": serialize({kind: {str(nu): preparation.initial_data(kind, 1000, nu) for nu in (1000, 3000, 4000)}
                                                       for kind in ("transverse", "longitudinal")}),
            "positive_normalized_preparation": "The complete frequency and slope changes are bounded by F/omega^5,G/omega^5, F=V_3+1,G=G_3+1. W_B>=omega/4. At u0=-1/2 use v=(2W_B)^(-1/2), v'=(-W_B'/(2W_B)-iW_B)v, then evolve the unchanged exact equation. The Cauchy Bogoliubov identities preserve CCR and give abs(B0)<=C_H/nu^6, C_H=10000000.",
            "Proca_reference_construction": "Use a_aux=1+chi(u+3)(a_actual-1), flat for u<=-2 and physical for u>=-1. It is smooth and >=1, hence globally hyperbolic on R times R³. A constrained positive Minkowski Proca seed propagates to this auxiliary spacetime, then from the shared Cauchy neighborhood to the actual one. The auxiliary geometry is not a physical model change.",
            "imported_theorem_and_hypotheses": "Moretti, Murro and Volpe, arXiv:2210.09278v3, Proposition 4.7: Proca time-slice extension and Hadamard propagation. Applied to the stated globally hyperbolic geometries, Cauchy neighborhoods and constant m²>0. https://arxiv.org/html/2210.09278v3 . No scalar state-of-low-energy theorem is imported.",
            "all_order_Hadamard_comparison": "On the compact auxiliary past interval, arbitrary fixed-order WKB references match the flat initial data and have mixing O(k^(-2N-1)). The new Cauchy data have the complete same WKB asymptotics. Their relative mixing with the propagated Proca Hadamard reference is faster than every inverse power. Physical vector reconstruction and its finite spacetime derivatives grow only polynomially on compact intervals, giving a smooth two-point difference. Positivity and CCR are retained, so the new state is Hadamard.",
            "physical_vector_reconstruction": "Transverse multiplier 1/sqrt(a); longitudinal spatial multiplier omega/(sqrt(a)*m) and temporal readout k(v'-d_L v)/(a^(5/2)*m*omega), up to Fourier phase. The constrained Minkowski polarization sum and positive three-mode Gram matrix are checked. The low-momentum new data equal the frozen data.",
            "quantitative_observable_transfer": "Uniformly on I and all k, each energy/pressure state-change integral is <=20C_H/m², each density derivative <=100C_H/m, and the clock-source difference <=100C_H/m+240C_H/m². Exact initial-state comparison supplies these constants; the qualitative auxiliary-reference constants are not inserted into the error budget.",
            "finite_state_change_regulator_limit": "Continue W_4(D) analytically on |D-3|<=1/4 and add the same D-independent Cauchy correction. For each analytic-sign transfer column, abs(A0-1)+abs(B0)<=C_H/nu^6. The continued combined readout difference is <=1000C_H*(25/16)/nu^5. Its dimensionally weighted radial tail is uniformly integrable, proving a holomorphic finite state-change limit and no new pole or finite counterterm. Compact-time distributional passage gives the clock-source limit.",
            "physical_scale_example": serialize(transfer.physical_bounds(10**12, 1000)),
            "physical_units_and_scope": "Restore tau^-4 for energy/pressure and tau^-5 for density derivatives and the physical phi=tau*u Euler source. Relative scales are M²/tau² and M²/tau³. At M*tau=10^12,m0*tau=1000 all five new-state totals remain below 10^-14. This is the stated vector-only Gaussian component at fixed background.",
            "controls": controls(),
            "verdict": "The explicit new vector state is Hadamard by constrained Proca comparison and has quantitatively transferred matched first-variation bounds. It does not solve the quantum-corrected light system or close original P8.",
            "not_established": ["Quantitative higher time/functional variations, other fields and interacting loops or omitted UV matching operators",
                                "Quantum-corrected light equations, constraints, principal cones, nonlinear stability or interacting cutoff",
                                "A common vacuum/UV parent, finite-gravity Regge remainder, V/G/B or original P8 closure"],
            "verification_boundary": "Exact symbolic and rational controls plus written all-order recurrence, cutoff, Proca Hadamard comparison and continuous transfer proofs using the explicitly identified propagation theorem; not proof-assistant formalization or a full quantum solution"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The controlled Hadamard vector-state report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.55.VECTOR_HADAMARD replay passed; new controlled vector state, original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
