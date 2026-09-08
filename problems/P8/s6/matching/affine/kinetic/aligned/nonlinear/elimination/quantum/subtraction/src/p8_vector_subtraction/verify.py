"""Read-only actual vector adiabatic subtraction and finite integral bound."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_vector_state import verify as parent

from . import controls as science_controls
from . import tail

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"finite-adiabatic-vector-integrals.json"
PARENT_SHA = "801b6d97564d15c23c4b5015be0ecbddaa8ee542b76cefce5b4aa04f33b157c8"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_vector_subtraction/*.py"))
            +sorted(ROOT.glob("tests/*.py"))+sorted(ROOT.glob("*.md"))
            +sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen physical vector energy comparison certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_50_fully_rebuilt": PARENT_SHA,
            "original_action_physical_frame_and_Gaussian_preparation_retained": True,
            "no_frozen_action_or_matched_finite_counterterm_changed": True}


def residuals():
    return {**tail.checks(), **science_controls.checks()}


def controls():
    bad = (True, False, sp.true, sp.false, 1.0, sp.Float(1), "1", sp.I,
           sp.oo, -sp.oo, sp.zoo, sp.nan, sp.Symbol("unproved"), 0, -1)
    calls = [lambda value=value: tail.physical_bounds(value, 1000) for value in bad]
    calls += [lambda value=value: tail.physical_bounds(10**12, value) for value in (*bad, 999)]
    calls += [lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan)]
    calls += [lambda value=value: science_controls.finite_moment((1-tail.wkb.z)**2, value)
              for value in (False, True, 0, 3, 1.0, sp.Integer(1), "1", sp.nan)]
    calls += [lambda value=value: science_controls.finite_moment(value, 2) for value in
              (sp.Integer(1), tail.wkb.z, sp.Float(1), sp.oo, sp.nan, 1.0)]
    calls += [lambda: tail.bounds("scalar")]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An invalid scale, mode, moment or report value was accepted")
    return {"rejected_inputs": rejected,
            "physical_clock_mass_derivatives_retained": True,
            "formal_derivative_order_not_naive_fixed_momentum_mass_order": True,
            "nonintegrable_moments_not_assigned_finite_Beta_integrals": True,
            "ordinary_Proca_benchmark_not_substituted_for_actual_lapse_energy": True,
            "fourth_order_subtraction_cannot_be_dropped": science_controls.incomplete_subtraction()["unsubtracted_fourth_order_energy_logarithmic_UV_coefficient_at_bounce"] != 0,
            "full_subtracted_integral_not_automatic_covariant_matching_or_VGB": True}


@cache
def build_report():
    pins = prior_checks()
    identities = residuals()
    exact = affine.certify_residuals(identities)
    proofs = science_controls.proof_checks()
    if not all(value is True for value in proofs.values()):
        raise ValueError("A continuous reference or integrated subtraction proof failed")
    generic = tail.algebra()
    return {"schema": 1, "claim": "P8-S6.51.VECTOR_SUBTRACTION", "date": "2026-09-08",
            "status": "FINITE_EXPLICITLY_SUBTRACTED_VECTOR_GAUSSIAN_ENERGY_AND_PRESSURE_CERTIFIED; ORIGINAL_P8_OPEN",
            "prior_sha256": pins,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md", "written_proofs": ["notes/proof.md"],
            "exact_residuals": exact, "named_exact_check_count": len(identities),
            "checked_scalar_entries": len(identities), "proof_checks": proofs,
            "literal_component": "Unchanged S6.42 retained vector, actual S6.50 physical lapse/pressure forms including a_N,b_N, and its exact normalized fourth-order Gaussian preparation on u in [-1/2,1/2], every real comoving momentum and m0*tau>=1000.",
            "physical_reference_weights": serialize(tail.physical_weights()),
            "derivative_order_prescription": "P2,P4,B2,B4,c1 scale with independent formal epsilon powers 2,4,3,5,1, respectively; omega,z,A,B have order zero. Subtract every term through epsilon^4 in the actual per-polarization omega*F/(4a_s³). This includes UV-finite parts of those full adiabatic coefficients.",
            "generic_reference_and_exact_tail": serialize({name: generic[name] for name in
                ("S", "shifted_rate", "actual", "zero", "second", "fourth", "remainder")}),
            "actual_full_subtraction_coefficients": serialize(tail.reference_terms()),
            "continuous_reference_tail_bounds": serialize({kind: tail.bounds(kind) for kind in ("transverse", "longitudinal")}),
            "reference_tail_and_momentum_integral": "abs(E(f)-E_AD[0:4])<=D/(4a_s³*omega^5) for each energy/pressure polarization, D=400000. integral d³k/(2pi)³ omega^-5=a_s³/(6pi²*m²). Summing all three polarizations gives <=D/(8pi²*m²)<D/(72m²), including the absolute integral.",
            "full_subtracted_exact_mode_integral": "Define the integral of E(v)-E_AD[0:4] by summing the absolutely integrable exact-mode/reference difference of S6.50 and the present reference tail. Uniformly on I, abs(rho_sub)<=12C/m+D/(72m²), abs(p_sub)<=108C/(5m)+D/(72m²), C=2000000, D=400000. This is the complete finite integral in this explicit subtraction prescription.",
            "physical_units_and_scheme_boundary": "Restore densities by tau^-4. Divide by M²/tau² for the displayed example. Arbitrary finite local terms in a different scheme are not controlled; no complete covariant mass/metric finite matching or all-order Hadamard claim is inferred.",
            "physical_scale_example": serialize(tail.physical_bounds(10**12, 1000)),
            "independent_ordinary_Proca_integrated_control": serialize(science_controls.ordinary_transverse()),
            "incomplete_subtraction_control": serialize(science_controls.incomplete_subtraction()),
            "incomplete_subtraction_consequence": "At u=0, 2F4_T+F4_L tends to -4/27 at high momentum. Omitting fourth order leaves -log(Lambda)/(54pi²) in the actual lapse energy. The finite exact evolution difference cannot cancel it.",
            "controls": controls(),
            "verdict": "Actual vector energy/pressure have a finite, explicitly prescribed all-momentum adiabatic subtraction with a quantified uniform bound for the selected Gaussian modes. Full covariant matching and original P8 remain open.",
            "not_established": ["Full covariant finite mass/metric counterterm matching, all-order Hadamard admissibility, all-time state control or a scheme-independent bound",
                                "Other field loops, higher-loop error, quantum-corrected bounce/constraints/cones, finite higher matching operators or interacting cutoff",
                                "Lorentz-invariant vacuum, finite-gravity Regge remainder, V/G/B admissibility or original P8 closure"],
            "verification_boundary": "Exact physical-weight comparisons, formal derivative-order series, exact rational tail, continuous coefficient envelopes, independent convergent radial moments and integrated Proca benchmark; written proof, not proof-assistant formalization"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The finite vector subtraction report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.51.VECTOR_SUBTRACTION replay passed; finite stated integral, original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
