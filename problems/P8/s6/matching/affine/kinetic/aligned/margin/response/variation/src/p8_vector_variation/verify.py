"""Read-only actual vector mode/readout variation and Cauchy-state boundary."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_margin_response import verify as parent

from . import modes, preparation, readouts

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"actual-vector-mode-variation.json"
PARENT_SHA = "a200c10d71ce6771c8ce921e120e0a00da13e53658ba75860d6e31d5737ae3b7"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_vector_variation/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen fixed-vector-stress response certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_57_fully_rebuilt": PARENT_SHA,
            "action_and_fixed_background_state_not_overwritten": True,
            "new_calculation_is_before_integrated_quantum_feedback_renormalization": True}


def residuals():
    return {**modes.checks(), **preparation.checks(), **readouts.checks()}


def controls():
    bad = (True, False, sp.true, sp.false, 1.0, sp.Float(1), "1", sp.I,
           sp.oo, -sp.oo, sp.zoo, sp.nan, sp.Symbol("unproved"), 0, -1, 1,
           sp.Rational(49, 100), sp.Rational(53, 100))
    calls = [lambda value=value: preparation.cauchy_mismatch(value) for value in bad]
    calls += [lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan)]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An invalid exact impedance diagnostic or report value was accepted")
    return {"rejected_inputs": rejected,
            "actual_temporal_mass_constraint_not_four_unconstrained_scalars": True,
            "original_canonical_variation_and_moving_oscillator_map_both_retained": True,
            "initial_covariance_term_not_set_to_zero_without_a_state_prescription": True,
            "energy_pressure_contact_terms_and_second_mass_jets_retained": True,
            "transverse_leading_conformal_control_not_full_Hadamard_claim": True,
            "off_clock_operator_not_silently_treated_as_ordinary_Proca": True,
            "fixed_Cauchy_prescription_failure_not_nonexistence_of_compatible_state": True}


@cache
def build_report():
    pins, identities = prior_checks(), residuals()
    exact, proofs = affine.certify_residuals(identities), preparation.proof_checks()
    if not all(value is True for value in proofs.values()):
        raise ValueError("An initial-force, impedance or causal-variation control failed")
    scalar_count = sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1 for value in identities.values())
    contacts = {kind: {name: data[name] for name in ("energy_contact_matrix", "pressure_contact_matrix")}
                for kind, data in readouts.matrices().items()}
    return {"schema": 1, "claim": "P8-S6.58.VECTOR_VARIATION", "date": "2026-09-08",
            "status": "ACTUAL_PER_MODE_CAUSAL_VARIATION_AND_MANDATORY_INITIAL_COVARIANCE_RESPONSE; ORIGINAL_P8_OPEN",
            "prior_sha256": pins,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md",
            "written_proofs": ["notes/modes.md", "notes/readouts.md", "notes/preparation.md"],
            "exact_residuals": exact, "named_exact_check_count": len(identities),
            "checked_scalar_entries": scalar_count, "proof_checks": proofs,
            "literal_component": "The unchanged retained vector on homogeneous physical histories with N>0,a>0 and the actual temporal/spatial mass functions. The temporal component is eliminated before canonical normalization. This derives per-mode response before UV integration, not a renormalized feedback norm.",
            "complete_mode_coefficients": serialize({kind: modes.canonical()[kind] for kind in ("transverse", "longitudinal")}),
            "actual_first_metric_variations": serialize(modes.first_variations()),
            "derivative_free_causal_phase_form": "For physical y=(w,p_w), M=[[0,g^-2],[-g² Omega_bare²,0]] and delta M=[[0,-2r/g²],[-g²(delta Omega_bare²+2r Omega_bare²),0]], r=delta log(g). Delta M depends algebraically on n and physical log-scale zeta. Delta y=Phi delta y_initial+integral Phi delta M y retains the initial-state term.",
            "equivalent_moving_oscillator_form": "The exact symplectic map C=[[g,0],[d*g,1/g]], d=g'/g, has delta C*C^-1=[[r,0],[r'+2d*r,-r]]. Its full variation gives delta Omega²=delta Omega_bare²-r''-2d*r'. The retarded kernel has the positive flat sine limit and the correct diagonal normalization. No artificial omission of the map or initial term is used.",
            "actual_readout_and_contact_definition": "At fixed physical canonical data and k_com, K_rho=(partial_N K_H)/a³ and K_p=-(a partial_a-2q partial_q)K_H/(3Na³). They reproduce the actual S6.50 clock readouts. Delta K=n partial_N K+zeta(a partial_a-2q partial_q)K includes changing volume and both second lapse mass jets.",
            "actual_mass_lapse_jets": serialize(readouts.actual_mass_jets()),
            "per_mode_local_contact_matrices": serialize(contacts),
            "covariance_and_observable_variation": "Delta Sigma'=M Delta Sigma+Delta Sigma M^T+delta M Sigma+Sigma delta M^T. Delta O=Tr(delta K Sigma+K Delta Sigma)/2. The initial covariance and the local contact terms are explicit. Their UV-divergent momentum integrals are not separately declared finite.",
            "actual_initial_force_bound": serialize(preparation.initial_force()),
            "nonzero_initial_geometry_change": "At L=10^12,R=1000 the exact local force polynomial exceeds 4R^4, pi²<10 and the full state error give F_n(-1/2)>5*10^-15. For the whole allowed margin range 0<J_e(-1/2)<1, so n_initial<-1/(4*10^14). The upper bound abs(n)<94*10^-14 keeps the original clock tube, with N<1 and p_affine>1/2.",
            "physical_Cauchy_impedance_test": "The transverse leading physical canonical impedance is k_com. The longitudinal impedance ratio on the exact initial chart lift has r_L²=a_m*b_m/(4p_affine²). Its exact rational factorization and a continuous polynomial sign bound prove 0<r_L²<1 for the actual p_affine in (1/2,21/40).",
            "initial_impedance_factorization_data": serialize(preparation.impedances()),
            "frozen_covariance_failure_scope": "Copying the old longitudinal physical coordinate/momentum covariance yields A_infinity=(r+1)/(2sqrt(r)), B_infinity=(r-1)/(2sqrt(r)), preserving CCR but with nonzero negative-frequency principal weight. Its vacuum-subtracted oscillator excitation energy has a quartically divergent radial tail. This fails the necessary leading adiabatic preparation condition; it does not exclude a compatible state family.",
            "distinct_oscillator_data_control": "Freezing normalized oscillator data is a different prescription. For transverse modes r_T=N exp(-omega), and r_T^4-1=(N²-1)[(h-1)N²+h]/h, so its leading frequency changes at every positive N!=1. The physical transverse conformal control does not establish all subleading state data.",
            "exact_impedance_examples": serialize({str(value): preparation.cauchy_mismatch(value)
                                                  for value in (sp.Rational(1, 2), sp.Rational(501, 1000))}),
            "renormalization_and_operator_boundary": "The contact and retarded state pieces must be combined with the corresponding full second-variation counterterms and a compatible initial covariance before a finite quantum feedback map is assigned. Off clock the anisotropic mass tensor changes the longitudinal principal operator; the ordinary-Proca propagation theorem is not silently transferred.",
            "controls": controls(),
            "verdict": "The actual per-mode vector and readout response is derived, and the actual fixed-source initial geometry is shown to require a nontrivial compatible covariance prescription. A renormalized integrated quantum response and original P8 remain open.",
            "not_established": ["A compatible all-order covariance family on the perturbed histories with quantitative error bounds",
                                "The complete renormalized retarded/contact second functional variation, feedback contraction or nonlinear quantum residual",
                                "A full coupled physical-cone, interacting-cutoff, common vacuum/UV parent, finite-gravity Regge, V/G/B or original P8 verdict"],
            "verification_boundary": "Exact symbolic/rational mode, contact, covariance and initial-force controls with written causal-ODE and continuous high-frequency preparation proofs; not proof-assistant formalization or a full quantum solution"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The actual vector variation report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.58.VECTOR_VARIATION replay passed; actual per-mode response, original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
