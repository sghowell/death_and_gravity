"""Read-only pinned replay of actual reference-stress and SEE-residual bounds."""

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp
from p8a_remainder import verify as prior

from . import bounds, independent, kernel, plateau, reconstruction

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"radiation-residual.json"
PRIOR_SHA = "cc4ec02bcb23e0ef01d0fb58ad4b95a3bb602a5bbbfc6871d044b2a09a78e963"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prior_checks():
    if sha(prior.REPORT) != PRIOR_SHA:
        raise ValueError("Pinned A.7 certificate changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return PRIOR_SHA


def verified_residuals(values):
    for name, value in values.items():
        if sp.simplify(value) != 0:
            raise ValueError(f"Nonzero actual-residual derivation identity: {name}")
    return dict.fromkeys(values, "0")


def verified_exclusion(name, value):
    value = sp.simplify(value)
    if value == 0:
        raise ValueError(f"Required actual-residual exclusion vanished: {name}")
    return str(value)


def checked_bounds():
    previous = json.loads(prior.REPORT.read_text())
    separate = independent.replay(previous)
    direct = bounds.calibration()
    for name in ("kernel_norm_constants", "effective_S_Born_norm_constants"):
        if list(map(str, direct[name])) != separate[name]:
            raise ValueError("Independent correct-clock kernel norm constants disagree")
    for name in ("Euler_first_constants", "derivative_first_constants", "R_squared_history_second_constants",
                 "first_constants", "second_constants"):
        if {key: str(value) for key, value in direct[name].items()} != separate[name]:
            raise ValueError(f"Independent actual-stress coefficients disagree: {name}")
    for name in ("rounded_first_constants", "rounded_second_constants"):
        if direct[name] != separate[name]:
            raise ValueError("Independent rounded coefficients disagree")
    delta = sp.Symbol("delta", nonnegative=True)
    records = {}
    for name in bounds.NAMES:
        first = direct["rounded_first_constants"][name]
        second = direct["rounded_second_constants"][name]
        stress_polynomial = delta*first+delta**2*second
        residual_polynomial = delta**2*(sp.Rational(bounds.FROZEN_NUMERATORS[name], 1922)
                                       + 2880*(first+delta*second))
        for polynomial in (stress_polynomial, residual_polynomial):
            if any(power[0] < 1 or coefficient <= 0 for power, coefficient in sp.Poly(polynomial, delta).terms()):
                raise ValueError("Finite-interval monotone positive majorant check failed")
        records[name] = {"actual_stress_error_polynomial": str(stress_polynomial),
                         "actual_SEE_residual_polynomial": str(residual_polynomial),
                         "positive_coefficients_and_monotonicity_on_closed_interval": True}
    example = separate["finite_physical_amplitude_example"]
    ratios = bounds.reference_accuracy(sp.Rational(example["delta_upper"]))
    if {name: str(value) for name, value in ratios.items()} != example["reference_stress_error_ratios"]:
        raise ValueError("Independent finite physical reference-scale ratios disagree")
    separate["monotone_amplitude_polynomials"] = records
    margins = kernel.elementary_margins()
    if any(value <= 0 for value in margins.values()):
        raise ValueError("A correct-clock local log cap was not proved")
    separate["correct_clock_kernel_rational_margins"] = {name: str(value) for name, value in margins.items()}
    plateau_margins = plateau.elementary_margins()
    for value in plateau_margins.values():
        if any(item <= 0 for item in (value if isinstance(value, list) else [value])):
            raise ValueError("An Euler or frozen-residual rational sign check failed")
    separate["plateau_rational_margins"] = {
        name: list(map(str, value)) if isinstance(value, list) else str(value)
        for name, value in plateau_margins.items()}
    return separate


def build_report():
    prior_hash = prior_checks()
    paths = sorted(ROOT.glob("src/p8a_residual/*.py"))+sorted(ROOT.glob("tests/*.py"))
    paths += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    eta = sp.Symbol("eta", positive=True)
    a, s, j, ell = (sp.Function(name)(eta) for name in ("a", "S", "J", "L"))
    hbar = sp.Symbol("hbar", positive=True)
    gamma, constant = sp.symbols("gamma C", real=True)
    exclusions = reconstruction.controls()
    preserved_wrong_constant = exclusions.pop("wrong_past_constant_can_still_be_conserved")
    return {
        "schema": 1, "claim": "P8-A.8", "date": "2026-09-06",
        "status": "FINITE_ACTUAL_REFERENCE_STRESS_AND_SEE_RESIDUAL_CERTIFIED; EXACT_SEE_QSEI_FOCUSING_OPEN",
        "prior_A7_sha256": prior_hash,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in paths},
        "formulation": "FORMULATION.md", "written_proof": "notes/proof.md", "source_boundary": "notes/sources.md",
        "scope": {
            "field": "one real free massless minimally coupled scalar",
            "state": "same positive zero-mean homogeneous isotropic quasifree Hadamard state transported from the unchanged radiation past",
            "same_A7_metric_family_and_exact_clock": True,
            "past_Cauchy_slice": "y=1/2, with active history starting at y=1 because all potential and remainder data vanish before preparation",
            "target": "2<=y<=3; compact interior targets retain the original near-target plateau",
            "FK_convention": "+---; R=+6*(Hdot+2H^2), GS Box equals FK Box; physical-sign residuals have opposite FK component signs",
            "named_prescription": "raw H_lambda with lambda=2*sqrt(2)*A*eta_star^2, physical gamma=0, no additional alpha*R Wick shift",
            "generic_symbolic_gamma_and_radiation_constant_retained": True,
            "generic_gamma_or_other_lambda_covered_by_numeric_constants": False,
            "absolute_stress_scheme_independent": False,
            "actual_past_data_fix_homogeneous_constant_to_zero": True,
            "actual_finite_amplitude_stress_vs_A3_reference_bound": True,
            "actual_finite_amplitude_SEE_residual_bound": True,
            "physical_epsilon_one_in_explicit_regime": True,
            "positive_components_for_specified_transported_state_in_accuracy_regime": True,
            "positive_EED_for_arbitrary_Hadamard_states": False,
            "source_model": "Einstein plus separately conserved ordinary radiation, zero Lambda, no additional explicit curvature-squared gravitational sources",
            "extra_curvature_sources_automatically_bounded": False,
            "exact_or_nearby_SEE_solution_or_stability": False,
            "self_consistent_past_initial_data_certified": False,
            "small_residual_on_entire_off_shell_preparation_interval": False,
            "second_order_corrected_metric_family_certified": False,
            "perturbed_QSEI_or_all_sampler_bound": False,
            "all_geodesic_focusing_or_P8a_completion": False,
        },
        "exact_trace_conservation_and_anomaly_identities": verified_residuals(reconstruction.identities()),
        "correct_clock_log_kernel_identities": verified_residuals(kernel.identities()),
        "exact_plateau_reference_and_frozen_residual_identities": verified_residuals(plateau.identities()),
        "generic_reconstructed_physical_stress": {name: str(value) for name, value in reconstruction.components(
            a, s, j, ell, eta, hbar=hbar, gamma=gamma, radiation_constant=constant).items()},
        "effective_anomaly_variable": "S=Rmode-K/2+U/10; W_eff=W-hbar*R/(240*pi^2); no scheme change",
        "state_fixed_histories": "J'=U'*S, L'=(a'/a)*U^2; both integrals and C/a^4 are fixed by the actual radiation past",
        "amplitude_domain": "0<=delta<=2048/190897521, delta=16*epsilon*d/(A^2*eta_star^4), d=kappa*hbar/(46080*pi^2)",
        "actual_stress_error_units": "hbar*(delta*C1+delta^2*C2)/(pi^2*A^4*eta_star^8)",
        "actual_SEE_residual": "delta^2/(A^2*eta_star^4)*[c_i/1922+2880*(C1_i+delta*C2_i)], c_i=(3,5,9)",
        "literal_1922_dictionary": "1922=2*31^2; 2048/(961*4096)=1/1922; not a power or typographical shorthand",
        "independently_replayed_finite_bounds": checked_bounds(),
        "negative_controls": {
            **{name: verified_exclusion(name, value) for name, value in exclusions.items()},
            "wrong_past_constant_can_still_be_conserved": preserved_wrong_constant,
            "required_initial_density_check_not_replaced_by_conservation": True,
            "input_rejections": "inexact/negative/nonfinite norms, missing jets, invalid amplitude and scales, unsupported coefficient selectors or alternate numeric scheme arguments",
        },
        "verification_boundary": [
            "A.7 and every earlier pinned state, prescription, all-frequency bound and source hash are replayed unchanged",
            "the analytic norm estimate differentiates the integrable retarded logarithm through the smooth flat-past source, not through a singular unregulated endpoint",
            "generic trace/conservation and finite tensor identities are symbolic; independent Laurent controls test nonzero-curvature algebra without claiming physical states",
            "a separate Fraction-only assembly from the pinned A.7 report verifies every numeric stress and residual coefficient and finite-interval ratio",
            "the certificate bounds an actual equation residual only on the specified target plateau; residual size is not an exact-solution or shadowing theorem",
            "positivity is proved only for the chosen transported state in the small-amplitude accuracy regime; full two-point QSEI and cosmological focusing remain open",
        ],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Radiation actual-residual certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8(a) A.8: actual reference-stress and finite SEE-residual bounds replay passed; exact SEE/QSEI/focusing OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
