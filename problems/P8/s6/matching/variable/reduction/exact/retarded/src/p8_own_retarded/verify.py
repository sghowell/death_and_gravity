"""Read-only finite-parameter causal own-f response certificate."""
import argparse
import hashlib
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

import sympy as sp
from p8_exact_stationary import verify as prior
from p8_variable_reduction.verify import serialize

from . import action, bounds, green

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"retarded-own-history.json"
PARENT_SHA = "727811a8563584cc511e749ff038eadaa416f02aae268a15fa4f90eb8ca4b6b3"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_files():
    return (sorted(ROOT.glob("src/p8_own_retarded/*.py"))
            +sorted(ROOT.glob("tests/*.py"))+sorted(ROOT.glob("*.md"))
            +sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(prior.REPORT) != PARENT_SHA:
        raise ValueError("The frozen exact stationary background report changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {"S6_33_recursively_rebuilt": PARENT_SHA,
            "same_unchanged_action_physical_metric_clock_and_free_chi": True,
            "uses_actual_finite_delta_own_f_background_not_a_limiting_substitute": True}


def interface_checks():
    coefficient, kernel = bounds.coefficient_box(), green.green_constants()
    values = {name: coefficient[name]-kernel[name] for name in ("k_min", "k_max", "s_min", "s_max")}
    values.update({"left": coefficient["x_min"]-kernel["left"],
                   "right": coefficient["x_max"]-kernel["right"],
                   "width": coefficient["width"]-kernel["width"],
                   "delta_max": coefficient["delta_max"]-green.response_constants()["delta_max"]})
    independently = action.independent_volterra_constants()
    values.update({"independent_R_x": independently["R_x"]-kernel["R_x_lower"],
                   "independent_R_ratio": independently["R_over_duration"]-kernel["R_over_distance_lower"],
                   "independent_Q0_floor": independently["Q0_lower_per_amplitude"]-green.response_constants()["Q_final_lower_over_eta"],
                   "independent_Qx0_floor": independently["Qx0_lower_per_amplitude"]-green.response_constants()["Q_x_final_lower_over_eta"]})
    return values


def numeric_response(delta=Fraction(1, 625), amplitude=1, parent_planck_squared=1, time_scale=1):
    """Strict int/Fraction calibration, not a numerical Green solver."""
    delta, eta, m2, tau = map(bounds.rational, (delta, amplitude, parent_planck_squared, time_scale))
    bounds.denominator_ratio(delta, 0)
    if eta <= 0 or m2 <= 0 or tau <= 0:
        raise ValueError("Positive amplitude, Planck-squared and time scale are required")
    response = green.response_constants()
    qlow = response["Q_final_lower_over_eta"]*eta
    qxlow = response["Q_x_final_lower_over_eta"]*eta
    root_delta = sp.sqrt(sp.Rational(delta.numerator, delta.denominator))
    c = 2+delta
    equation_floor = m2*128*qlow/(tau**2*delta*c)
    return {"delta": delta, "c": c, "amplitude": eta, "M_squared": m2, "tau": tau,
            "Q0_lower": qlow, "Qx0_lower": qxlow,
            "QT0_lower": qxlow/(tau*root_delta),
            "physical_window_duration": tau*root_delta/4,
            "germ_local_minimax_Q0_error_lower": qlow/2,
            "normalized_original_g_equation_magnitude_lower": 128*qlow/c,
            "physical_original_g_equation_magnitude_lower": equation_floor,
            "literal_g_action_derivative_magnitude_lower": equation_floor/2,
            "E_g_convention": "-2 times the literal norm-two action Euler derivative",
            "fixed_physical_low_frequency_band_claim": False,
            "canonical_matter_source_or_undriven_physical_g_solution": False}


def controls():
    bad = (True, False, 1.0, sp.Float("1"), sp.oo, -sp.oo, sp.nan, "1/625", sp.Symbol("delta"))
    calls = [lambda value=value: numeric_response(delta=value) for value in bad]
    calls += [lambda value=value: numeric_response(delta=value) for value in (0, -1, Fraction(1, 500))]
    calls += [lambda: numeric_response(amplitude=0), lambda: numeric_response(parent_planck_squared=0),
              lambda: numeric_response(time_scale=0)]
    calls += [lambda value=value: green.volterra_bounds(width=value) for value in (True, 0.25, "1/4", sp.oo, sp.nan)]
    calls += [lambda params=params: green.volterra_bounds(**params)
              for params in ({"k_min": 0}, {"k_min": 5, "k_max": 4}, {"s_max": -1}, {"width": -1})]
    calls += [lambda value=value: green.flat_derivative_polynomials(value) for value in (True, 1.0, "1", sp.Float(1), -1)]
    calls += [lambda value=value: serialize(value) for value in (0.1, sp.Float("0.1"), sp.oo, sp.nan, sp.I)]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An invalid physical/Green/arithmetic input was accepted")
    if green.volterra_bounds(width=Fraction(1, 2))["positive_flux_bootstrap"]:
        raise ValueError("A longer unproved positive-kernel interval was promoted")
    return {"rejected_inputs": rejected, "longer_interval_bootstrap_rejected": True,
            "zero_hidden_data_at_left_endpoint_explicit": True,
            "flux_jump_is_one_not_R_derivative_one": True,
            "whole_variable_kinetic_derivative_retained": True,
            "prescribed_metric_is_not_a_conserved_matter_source": True,
            "shrinking_duration_is_not_a_fixed_low_frequency_class": True,
            "causal_equation_not_silently_a_single_copy_variational_action": True,
            "literal_action_derivative_factor_two_recorded": True}


@cache
def build_report():
    pins = prior_checks()
    for name in ("test_own_retarded_action_audit.py", "test_own_retarded_bounds.py", "test_own_retarded_green.py"):
        if not (ROOT/"tests"/name).is_file():
            raise ValueError("Each separately authored scientific audit must be pinned")
    groups = {"literal_time_action_and_original_metric_readout": action.identities(),
              "exact_Green_and_pulse_constants": green.identities(),
              "independent_interfaces_and_Green_recalculation": interface_checks()}
    if any(value != 0 for group in groups.values() for value in group.values()):
        raise ValueError("A literal action, Green constant or interface identity failed")
    proof_checks = {"coefficient_box": bounds.checks(), "Green_and_pulse": green.checks()}
    if not all(value is True for group in proof_checks.values() for value in group.values()):
        raise ValueError("A continuous coefficient or positive-Green proof check failed")
    return {
        "schema": 1, "claim": "P8-S6.34.OWN_RETARDED", "date": "2026-09-07",
        "status": "FINITE_DELTA_RETARDED_OWN_F_HISTORY_AND_ORIGINAL_G_EQUATION_READOUT; ORIGINAL_P8_OPEN",
        "prior_sha256": pins,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
        "formulation": "FORMULATION.md", "written_proofs": ["notes/bounds.md", "notes/green.md"],
        "source_audit": "notes/sources.md",
        "domain": "0<delta<=1/625; x=T/(tau sqrt(delta)) in[-1/4,0]; M,tau>0",
        "input_state_contract": "zero history OR prescribed smooth0<=q<=eta supported inside(-1/4,-1/8), integral(q)>=eta/16; Q=Q_x=0 at-1/4; positive lower bounds concern the pulse subclass",
        "exact_residuals": {name: dict.fromkeys(group, "0") for name, group in groups.items()},
        "exact_residual_count": sum(map(len, groups.values())),
        "proof_checks": proof_checks,
        "finite_delta_coefficient_box": serialize(bounds.calibration()),
        "Green_and_pulse_calibration": serialize(green.calibration()),
        "literal_original_g_readout": serialize(action.readout()),
        "independent_Green_inequalities": serialize(action.independent_volterra_constants()),
        "physical_numeric_calibration": serialize(numeric_response(amplitude=Fraction(1, 100), parent_planck_squared=3, time_scale=2)),
        "germ_local_minimax_Q0_error_per_amplitude": "63/2048",
        "controls": controls(),
        "not_established": [
            "An undriven physical-g/clock or full nonlinear parent solution",
            "A prescribed conserved matter stress, free-chi source loading or full source matching",
            "A fixed physical temporal band, cutoff, or general controlled-EFT exclusion",
            "A nonlinear hidden response, full covariant CD coefficient match or scalar health",
            "A variational single-copy action obtained merely by inserting a retarded inverse",
            "Positivity/V/G applicability, UV completion or original P8 closure"],
        "verification_boundary": "Written continuous finite-parameter coefficient and first-zero/Volterra proof, exact margins and three separately authored scientific audits; not proof-assistant formalization",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The retarded own-f history report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.34.OWN_RETARDED replay passed; original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
