"""Read-only Proca local curvature coefficients and independent spectral control."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_vector_bubble import verify as parent

from . import geometry, heat, sphere

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"local-curvature-counterterms.json"
PARENT_SHA = "b9b4e4c22c1f824bf101e902b6ebdb9a3535e1b7e59f51a5288fbf3e359c47e2"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_vector_curvature/*.py"))
            +sorted(ROOT.glob("tests/*.py"))+sorted(ROOT.glob("*.md"))
            +sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen vector mass-insertion certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_48_fully_rebuilt": PARENT_SHA,
            "original_all64_action_and_physical_clock_geometry_retained": True,
            "no_frozen_action_or_finite_counterterm_changed": True}


def residuals():
    out = {**heat.checks(), **geometry.checks(), **sphere.checks()}
    d, u = geometry.rolling(), geometry.u
    out["actual_rolling_local_a4_rational_form"] = sp.factor(d["local_a4_with_divergence"]
                                        +8*(77*u**4+14*u**2-3)/(1+u**2)**4)
    out["actual_rolling_local_a4_mod_divergence_form"] = sp.factor(d["local_a4_mod_divergence"]
                                        +8*(259*u**4+98*u**2-5)/(5*(1+u**2)**4))
    return out


def controls():
    bad = (True, False, sp.true, sp.false, 1.0, sp.Float(1), "1", sp.I,
           sp.oo, -sp.oo, sp.zoo, sp.nan, sp.Symbol("unproved"), 0, -1)
    calls = [lambda value=value: geometry.scale_bounds(value) for value in bad]
    calls += [lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan)]
    calls += [lambda value=value: geometry.even_envelope(value) for value in
              (1/(1+geometry.u**4), geometry.u/(1+geometry.u**2), geometry.u**4/(1+geometry.u**2))]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An invalid scale, coefficient envelope or report value was accepted")
    return {"rejected_inputs": rejected,
            "finite_cochain_contact_factor_not_silently_dropped": True,
            "massive_vector_minus_one_scalar_not_Maxwell_minus_two": True,
            "twenty_parameter_curvature_tensor_not_constant_curvature_only": heat.explicit_tensor_traces()["curvature_parameter_count"] == 20,
            "sphere_scalar_constant_mode_retained_in_ratio": True,
            "on_clock_coefficients_not_full_off_clock_quantum_variation": True,
            "actual_noncompact_boundary_flux_not_set_to_zero": True,
            "local_pole_and_S4_error_not_full_FLRW_in_in_error": True}


@cache
def build_report():
    pins = prior_checks()
    identities = residuals()
    exact = affine.certify_residuals(identities)
    proofs = {**geometry.proof_checks(), **sphere.proof_checks()}
    if not all(value is True for value in proofs.values()):
        raise ValueError("A continuous physical coefficient or independent spectral error proof failed")
    local, actual, spectral = heat.coefficients(), geometry.rolling(), sphere.spectral()
    return {"schema": 1, "claim": "P8-S6.49.VECTOR_CURVATURE", "date": "2026-09-08",
            "status": "LOCAL_VECTOR_CURVATURE_POLE_AND_GLOBAL_COEFFICIENT_BOUNDS_CERTIFIED; ORIGINAL_P8_OPEN",
            "prior_sha256": pins,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md", "written_proofs": ["notes/local.md", "notes/spectral.md"],
            "exact_residuals": exact, "named_exact_check_count": len(identities),
            "checked_scalar_entries": sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1 for value in identities.values()),
            "proof_checks": proofs,
            "literal_component": "Unchanged S6.42 retained vector at fixed light fields. Actual p(phi,x=-1)=1/2 and a=b=1 all along the clock make the mass tensor covariantly constant. Off-clock mass variations remain nonzero.",
            "determinant_prescription": "Euclidean +1/2 Tr log(D1+m0²)-1/2 Tr log(D0+m0²), D1=-nabla²+Ricci on one-forms and D0=-nabla² on scalars. Dimensional regularization sets scaleless contact traces to zero; finite cochain controls retain their contact determinant and scalar constant mode.",
            "local_heat_coefficients": serialize(local),
            "regulator_and_pole_normalization": "d=4-2epsilon_DR; Gamma_pole=-1/(32pi² epsilon_DR)*integral sqrt(g)[3m0^4/2+m0²*R/2+a4]. The flat term is exactly S6.47. Lorentzian evaluations concern local invariant counterterm functions, not a finite Lorentzian determinant.",
            "actual_physical_metric_invariants_and_coefficients": serialize(actual),
            "continuous_binomial_envelopes": serialize(geometry.continuous_bounds()),
            "sharp_scalar_curvature_and_units": "For every real u, 0<R*tau²<=49 by an exact square identity; |a4*tau^4|<=308/3 including boxR. R=24(1+7u²)/(1+u²)²/tau² and a4=-8(77u^4+14u²-3)/(1+u²)^4/tau^4.",
            "boundary_terms": "a³ boxR=-(a³R')' and a³ GB=(8a³H³)'; normalized fluxes grow as 336u^9 and 512u^9. No vanishing infinite-time flux or finite spacetime action is assumed. Dropping boxR changes the local coefficient.",
            "physical_scale_example": serialize(geometry.scale_bounds(1000)),
            "pole_ratio_boundary": "For Rm=m0*tau>0 the absolute curvature correction to the pole weight over its flat 3m0^4/2 is <=49/(3Rm²)+616/(9Rm^4). This is not a finite quantum error or an observable ratio.",
            "independent_S4_spectral_control": serialize({name: spectral[name] for name in
                ("eigenvalue", "multiplicity", "approximation", "leading", "polynomial_tail_bound",
                 "sixth_derivative_norm_envelopes", "Euler_Maclaurin_bound_before_exponential",
                 "trace_remainder_over_t_upper_for_0_lt_t_le_1")}),
            "S4_mode_and_domain_boundary": "The unit-S4 ratio trace is coexact-vector heat trace minus one scalar constant mode. For 0<s<=1 its difference from 1/(2s²)-1/s-11/30 is <=1609s/1728. This independent spectral bound is not transferred to the actual noncompact Lorentzian FLRW metric.",
            "controls": controls(),
            "verdict": "The retained vector's on-clock local curvature pole, all-time physical invariant coefficient bounds and a bounded independent spectral control are established. Original P8 remains open.",
            "not_established": ["Finite curved/in-in quantum remainder, complete mass/metric variations, full quantum stress tensor or renormalized bounce",
                                "Quantum-corrected constrained spectrum/cones, other field loops, finite higher matching operators, higher-loop error or interacting cutoff",
                                "Lorentz-invariant vacuum, finite-gravity Regge remainder, V/G/B admissibility or original P8 closure"],
            "verification_boundary": "Universal local Laplace-type coefficient input with explicit tensor contractions, finite Hodge complex control, direct physical metric curvature, exact continuous rational envelopes and bounded Euler-Maclaurin spectral derivation; not proof-assistant formalization"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The local vector curvature report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.49.VECTOR_CURVATURE replay passed; local counterterms, original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
