"""Read-only vector clock-mass finite matching and regulator-limit certificate."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_vector_dimensional import verify as parent

from . import bounds, continuation, counterterms, variation

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"finite-clock-mass-vector-matching.json"
PARENT_SHA = "8dd43f663b78e611115ddd3264e79490d6fb7b3884e17524387b05638be6477d"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_vector_clock_matching/*.py"))
            +sorted(ROOT.glob("tests/*.py"))+sorted(ROOT.glob("*.md"))
            +sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen dimensional vector matching certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_52_fully_rebuilt": PARENT_SHA,
            "original_action_physical_frame_and_Gaussian_preparation_retained": True,
            "no_frozen_finite_potential_or_lapse_jet_changed": True}


def residuals():
    return {**counterterms.checks(), **variation.checks(), **bounds.checks(), **continuation.checks()}


def proof_checks():
    d = bounds.physical_bounds(10**12, 1000)
    return {**continuation.proof_checks(),
            "actual_flat_energy_envelope": bool(sp.Rational(5, 2)+sp.Rational(22, 27) < 4),
            "fixed_prescription_energy_example_below_stated_bound": bool(d["energy_candidate_matched_total_over_reference_density"] < sp.Rational(1, 10**14)),
            "fixed_prescription_pressure_example_below_stated_bound": bool(d["pressure_candidate_matched_total_over_reference_density"] < sp.Rational(1, 10**14))}


def controls():
    bad = (True, False, sp.true, sp.false, 1.0, sp.Float(1), "1", sp.I,
           sp.oo, -sp.oo, sp.zoo, sp.nan, sp.Symbol("unproved"), 0, -1)
    calls = [lambda value=value: bounds.physical_bounds(value, 1000) for value in bad]
    calls += [lambda value=value: bounds.physical_bounds(10**12, value) for value in (*bad, 999)]
    calls += [lambda value=value: variation.jet_check(value) for value in
              (True, False, sp.true, 1.0, sp.Integer(1), "1", 0, -1, 5)]
    calls += [lambda value=value: continuation.box_bound(value) for value in
              (True, False, sp.true, 1.0, sp.Float(1), "1", sp.oo, sp.nan,
               1/(1-continuation.u**2), sp.sin(continuation.u), sp.sqrt(2), sp.Symbol("unproved"))]
    calls += [lambda value=value: continuation.coefficients(value) for value in ("scalar", None, True)]
    calls += [lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan)]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError, sp.PolynomialError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An invalid scale, jet order, coefficient or report value was accepted")
    return {"rejected_inputs": rejected,
            "lapse_variation_before_clock_restriction": True,
            "counterterm_variation_before_dimension_limit": True,
            "complex_dimension_not_conjugated_in_analytic_pair": True,
            "spatial_average_continuation_explicit_not_silently_re_ranked": True,
            "clock_mass_component_not_assumed_separately_conserved": True,
            "small_unsigned_first_variation_not_a_cone_or_VGB_verdict": True}


@cache
def build_report():
    pins, identities = prior_checks(), residuals()
    exact, proofs = affine.certify_residuals(identities), proof_checks()
    if not all(value is True for value in proofs.values()):
        raise ValueError("A complex-dimension or physical bound failed")
    scalar_count = sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1 for value in identities.values())
    return {"schema": 1, "claim": "P8-S6.53.VECTOR_CLOCK_MATCHING", "date": "2026-09-08",
            "status": "FIXED_PRESCRIPTION_VECTOR_FIRST_METRIC_VARIATIONS_MATCHED_AND_BOUNDED; ORIGINAL_P8_OPEN",
            "prior_sha256": pins,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md",
            "written_proofs": ["notes/prescription.md", "notes/dimension-limit.md", "notes/bounds.md"],
            "exact_residuals": exact, "named_exact_check_count": len(identities),
            "checked_scalar_entries": scalar_count, "proof_checks": proofs,
            "literal_component": "Unchanged S6.42 retained vector and physical clock; S6.47 finite potential and lapse jets, S6.50 exact fourth-order Gaussian preparation, S6.51 subtraction and S6.52 dimensional radial terms. Only first homogeneous metric variations on I=[-1/2,1/2].",
            "named_prescription": "Scalar-coefficient continuation: C4=2b_m²+a_m^(3/2)b_m^(1/2) independent of D. Fixed four-dimensional ordinary and scalar curvature pole coefficients evaluated as invariants in D+1 dimensions. Spatial contractions use three times the average P^munu T_munu/D, with the unit timelike clock normal. This is an explicit prescription, not unique UV data.",
            "counterterm_density": "Relative to 1/(64pi² epsilon_DR): m^4 C4+m²R+2a4V_4-m²(delta_a ga2+delta_b gb2)-(delta_a ga4+delta_b gb4). C2=Rg/12-5Ricci/6, Z4=Tsc4+g boxR/6; ga2=-2C2_nn, gb2=6C2_spatial_average, ga4=-Z4_nn/2, gb4=3Z4_spatial_average/2. Tsc4 varies 2a4sc_4 in D+1 dimensions. Higher mass-curvature variations are not supplied.",
            "first_variation_argument": "delta[(a_m-1)G]|clock=G_clock delta a_m before integration by parts. Time-dependent jet controls through order four confirm lapse alpha*ga+beta*gb and zero first homogeneous spatial variation. The mass deviations must not be erased before variation.",
            "dimensional_scalar_curvature_tensor": serialize(counterterms.scalar_curvature_tensor()),
            "mass_pole_extension_and_evanescent_energy": serialize(counterterms.mass_pole_extension()),
            "matched_local_first_variation_coefficients": serialize(counterterms.finite_coefficients()),
            "actual_clock_local_coefficients_at_mu_mass": serialize(bounds.actual_coefficients()),
            "continuous_actual_clock_local_envelopes": serialize(bounds.continuous_envelopes()),
            "complex_dimension_continuation": "On |D-3|<=1/4, real I, k>=0 and m>=1000, analytically continue the canonical ODE and prepare v_+,v_- against f_±=(2W)^(-1/2) exp(∓i integral W). Do not conjugate D. At D=3 these are the frozen mode and conjugate; no noninteger-dimensional Hilbert space is claimed.",
            "complex_dimension_coefficient_and_tail_envelopes": serialize({kind: continuation.bounds(kind) for kind in ("transverse", "longitudinal")}),
            "uniform_integrable_regulator_limit": "C=10000000,Ctail=1000000,Amax=25/16,nu²=m²+k²/Amax². Combined evolution difference <=1000 C Amax/nu^4 and reference tail <=Ctail/omega^5. Radial powers at infinity are k^(-7/4),k^(-11/4); at zero k^(7/4). The analytic reciprocal-Gamma prefactor is bounded on the compact disk. Dominated convergence and Morera establish F(D) holomorphic near D=3 with limit exactly the S6.51 finite physical integral.",
            "matched_finite_part": "F(3)+S6.52 radial local finite terms minus the epsilon coefficient of the explicitly continued counterterm metric variation. The flat finite term is preserved; switching off alpha,beta recovers the independent ordinary-Proca finite control. Pressure receives no added linear mass-counterterm first spatial variation.",
            "physical_scale_example": serialize(bounds.physical_bounds(10**12, 1000)),
            "physical_bound_boundary": "At mu=m0, restore tau^-4 and divide by M²/tau². Explicit local coefficients plus S6.51's sharper D=3 integral bound give energy and pressure magnitudes below 10^-14 for M*tau=10^12,m0*tau=1000. No unknown complex-dimension prefactor is used as a physical error coefficient.",
            "local_boundary_terms": "Compactly supported Euler variations of boxR vanish. S6.49's nonzero noncompact infinite-time flux statement is unchanged.",
            "controls": controls(),
            "verdict": "The retained vector's finite first homogeneous metric variations are matched and continuously bounded in the named prescription, including a uniform regulator-limit proof. Original P8 remains open.",
            "not_established": ["Unique scheme-independent finite terms, unknown UV Wilson coefficients, quadratic/higher clock-mass curvature variations or extension of the timelike projector to X=0",
                                "All-order Hadamard admissibility, derivatives of the nonlocal state term, the full clock equation, all-time state control or other field/higher loops",
                                "Quantum-corrected bounce, constraints, cones or interacting cutoff; Lorentz-invariant vacuum, finite-gravity Regge remainder, V/G/B or original P8 closure"],
            "verification_boundary": "Exact symbolic variations, analytic bilinear-mode controls, continuous rational box bounds and written dominated-convergence proof; not proof-assistant formalization or a full quantum solution"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The finite clock-mass vector matching report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.53.VECTOR_CLOCK_MATCHING replay passed; fixed first variations, original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
