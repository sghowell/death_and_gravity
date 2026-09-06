"""Read-only certificate replay for the actual all-H2-sampler QSEI."""

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp
from p8a_residual import bounds as stress_bounds
from p8a_residual import verify as prior

from . import independent, mode_bounds, sampling

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"prepared-qsei.json"
PRIOR_SHA = "dbd88193542224b3eeb76e18923744f3317791ad46961caca3562a1525ba56d1"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prior_checks():
    if sha(prior.REPORT) != PRIOR_SHA:
        raise ValueError("The pinned A.8 actual-stress certificate changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return PRIOR_SHA


def zero_checks(values):
    for name, value in values.items():
        if sp.simplify(value) != 0:
            raise ValueError(f"Nonzero all-sampler identity: {name}")
    return dict.fromkeys(values, "0")


def serialized(values):
    if isinstance(values, dict):
        return {key: serialized(value) for key, value in values.items()}
    if isinstance(values, list):
        return list(map(serialized, values))
    return str(values)


def checked_constants():
    modes, sampler = mode_bounds.calibration(), sampling.calibration()
    separate = independent.replay(json.loads(prior.prior.REPORT.read_text()))
    for key, value in {**modes, **sampler}.items():
        if key in separate and serialized(value) != separate[key]:
            raise ValueError(f"Fraction-only QSEI replay disagrees: {key}")
    accuracy = stress_bounds.reference_accuracy(sampling.DELTA_MAX)
    if any(value >= sp.Rational(1, 100) for value in accuracy.values()):
        raise ValueError("The absolute bound lacks its positive actual reference EED")
    separate["A8_actual_reference_EED_error_ratio"] = str(accuracy["EED"])
    # Independent rational geometric margins on s in [2,9/2], delta<=1/2.
    margins = {"inverse_a_squared": sp.Rational(31, 2)-4,
               "Hubble_cap_identity": sp.Rational(8, 31)-1/(4*sp.Rational(31, 32)),
               "Hubble_derivative_cap_identity": sp.Rational(132, 961)
                   -sp.Rational(33, 32)/(8*sp.Rational(31, 32)**2),
               "rational_sqrt_two_upper_margin": sp.Rational(3, 2)**2-2,
               "spectral_cross_triangle_margin": 8-sp.Rational(8, 3)**2}
    if margins["inverse_a_squared"] <= 0 or margins["rational_sqrt_two_upper_margin"] <= 0:
        raise ValueError("A rational geometry or square-root margin failed")
    if margins["spectral_cross_triangle_margin"] <= 0:
        raise ValueError("The auxiliary norm triangle coefficient is too small")
    if any(margins[name] != 0 for name in ("Hubble_cap_identity", "Hubble_derivative_cap_identity")):
        raise ValueError("A plateau Hubble constant is inconsistent")
    separate["rational_geometry_and_norm_margins"] = serialized(margins)
    return {"phase_resolved_modes": serialized(modes), "all_sampler_constants": serialized(sampler),
            "independent_Fraction_replay": separate}


def build_report():
    previous = prior_checks()
    paths = sorted(ROOT.glob("src/p8a_qsei/*.py"))+sorted(ROOT.glob("tests/*.py"))
    paths += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {"schema": 1, "claim": "P8-A.9", "date": "2026-09-06",
            "status": "ACTUAL_PREPARED_METRIC_ALL_SAMPLER_QSEI_CERTIFIED; SEE_AND_COSMOLOGICAL_FOCUSING_OPEN",
            "prior_A8_sha256": previous,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in paths},
            "formulation": "FORMULATION.md", "written_proof": "notes/proof.md",
            "source_boundary": "notes/sources.md",
            "scope": {"same_actual_A8_metric_state_preparation_and_named_scheme": True,
                      "field": "one free real massless minimally coupled scalar",
                      "reference": "actual transported positive homogeneous isotropic Hadamard state; not a Born state",
                      "target_states": "all Hadamard target states; no target homogeneity, quasifree, mean or energy restriction",
                      "samplers": "all real smooth compact samplers in 2<y<3, extending by density to H2_0 on each compact interval",
                      "delta_domain": "0<=delta<=10^-14",
                      "proper_time_bound": "integral h^2 E_omega dt >= -5*hbar/(16*pi^2)*integral |h_ddot|^2 dt",
                      "coefficient_is_sharp": False,
                      "finite_reference_credit_uses_A8_named_raw_H_lambda_gamma_zero": True,
                      "alternative_finite_prescription_automatically_covered": False,
                      "difference_QSEI_renormalization_independent": True,
                      "auxiliary_flat_kernel_is_a_perturbed_state": False,
                      "IR_and_UV_and_both_Fourier_frequencies_controlled": True,
                      "third_sampler_derivative_required": False,
                      "all_comoving_lines_covered_on_target": True,
                      "SEE_solution_or_nearby_metric_control": False,
                      "complete_A1_geodesic_domain_or_initial_curvature_hypothesis": False,
                      "massive_interacting_or_all_realistic_fields": False,
                      "new_cosmological_incompleteness_or_full_P8a_closure": False},
            "exact_phase_mode_identities": zero_checks(mode_bounds.identities()),
            "exact_spectral_moments": zero_checks(sampling.spectral_identities()),
            "exact_actual_proper_clock_identities": zero_checks(sampling.clock_identities()),
            "derived_and_independent_constants": checked_constants(),
            "negative_controls": {"nonflat_past_derivative_transfer": "boundary exp(2ik(x-1))*(U*b)^(j)(1) need not vanish",
                                  "wrong_complex_half_Parseval": "half of the full Fourier norm is not an upper bound for a general complex mode-error product",
                                  "UV_zeroth_order_absolute_mode_integral": "loses positive-frequency decay and is divergent",
                                  "auxiliary_modes_as_actual_state": "KG residual U*exp(-ikx) is nonzero",
                                  "unbounded_reference_credit": "difference bound alone is not the displayed absolute inequality",
                                  "input_rejections": "inexact, negative, nonfinite or above-cap amplitudes/norms; no arbitrary scheme or metric supplied to fixed-model API"},
            "verification_boundary": ["The analytic positive-type argument and Volterra/Sobolev proof are written, not Lean formalized",
                                      "All numerical mode and sampler constants are exact rationals with independent polynomial-Fraction replay",
                                      "Both inverse-frequency and infrared estimates are for actual modes, not only a coincident stress remainder",
                                      "The prepared metric is still off shell; neither the original full P8(a) statement nor P8 overall is closed"]}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The prepared all-sampler QSEI certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8(a) A.9: actual all-H2-sampler QSEI replay passed; SEE and cosmological focusing OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
