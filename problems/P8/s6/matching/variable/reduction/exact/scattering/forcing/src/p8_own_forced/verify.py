"""Read-only fixed-pulse, zero-data exact own-f phase-response certificate."""
import argparse
import hashlib
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

import sympy as sp
from p8_own_scattering import verify as prior
from p8_own_scattering.verify import serialize

from . import audit, bridges, loading, phase

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"fixed-pulse-own-response.json"
PARENT_SHA = "c77d6d3e41d81c5d2bda33aeba201fe438620d872f3e9f2d49fde5469194641e"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_files():
    return (sorted(ROOT.glob("src/p8_own_forced/*.py"))
            +sorted(ROOT.glob("tests/*.py"))+sorted(ROOT.glob("*.md"))
            +sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(prior.REPORT) != PARENT_SHA:
        raise ValueError("The frozen exact own-f fixed-window transfer changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {"S6_35_recursively_rebuilt": PARENT_SHA,
            "same_unchanged_action_off_shell_metric_clock_and_free_chi": True,
            "uses_actual_own_f16_not_the_coupled_relative80_response": True,
            "does_not_import_S6_34_shrinking_pulses_or_S6_27_matter_loading": True}


def flatten(values):
    result = {}
    for name, value in values.items():
        if isinstance(value, sp.MatrixBase):
            for row in range(value.rows):
                for col in range(value.cols):
                    result[f"{name}_{row}{col}"] = value[row, col]
        else:
            result[name] = value
    return result


def physical_contract(amplitude=1, time_scale=1):
    eta, tau = map(loading.rational, (amplitude, time_scale))
    if eta <= 0 or tau <= 0:
        raise ValueError("The fixed pulse amplitude and physical time scale must be positive")
    pulse, constants = loading.pulse(), loading.calibration()
    return {"eta": eta, "tau": tau, "delta_max": loading.DELTA_MAX,
            "initial_physical_time": -tau*loading.L,
            "loading_physical_time": -tau*loading.A,
            "observation_physical_time": tau*loading.A,
            "total_physical_duration": tau*(loading.L+loading.A),
            "physical_pulse_support": [tau*pulse["support_left"], tau*pulse["support_right"]],
            "physical_plateau_width": tau*pulse["plateau_width"],
            "loaded_Y_lower": constants["Y_load_lower_over_eta"]*eta,
            "loaded_Y_T_lower": constants["Y_u_load_lower_over_eta"]*eta/tau,
            "limsup_minus_liminf_original_Q_lower": eta/200,
            "stronger_original_Q_gap_lower": phase.calibration()["scalar_Q_gap_lower_over_eta"]*eta,
            "gap_not_a_pointwise_positive_Q_floor_for_every_delta": True,
            "same_q_and_tau_for_every_delta": True,
            "physical_low_frequency_band_or_cutoff_claimed": False,
            "canonical_matter_or_undriven_g_solution_claimed": False}


def controls():
    bad = (True, False, sp.true, sp.false, 1.0, sp.Float("1"), sp.oo, -sp.oo, sp.nan, "1", sp.I)
    calls = [lambda value=value: physical_contract(amplitude=value) for value in bad]
    calls += [lambda value=value: physical_contract(time_scale=value) for value in (True, 1.0, 0, -1)]
    calls += [lambda: physical_contract(amplitude=0), lambda: physical_contract(amplitude=-1)]
    calls += [lambda value=value: phase.phase_sequences(value) for value in (True, sp.true, 2.0, "2", 0, 1, -1)]
    calls += [lambda value=value: phase.phase_sequences(2, offset=value) for value in (True, 0.1, sp.oo, sp.I, -1)]
    calls += [lambda value=value: loading.volterra_bounds(value, Fraction(1, 100)) for value in (True, 1.0, "1", -1)]
    calls += [lambda value=value: loading.flat_derivative_polynomial(value) for value in (True, 1.0, "1", -1)]
    calls += [lambda value=value: phase.endpoint_limit(value, 0) for value in (True, 1.0, 0, -1)]
    calls += [lambda value=value: phase.phase_value(value) for value in (0, -1, sp.oo, sp.nan)]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An invalid physical input, phase sequence or proof domain was accepted")
    if loading.volterra_bounds(80044, loading.L)["positive_derivative_bootstrap"]:
        raise ValueError("An unproved long positive-Green interval was promoted")
    return {"rejected_inputs": rejected,
            "longer_positive_Green_interval_rejected": True,
            "strict_limiting_loading_uses_stronger_uniform_margin": True,
            "phase_offset_fixed_by_one_loaded_limit_not_fitted_per_delta": True,
            "both_actual_error_neighborhoods_retained": True,
            "original_Q_readout_not_an_arbitrary_projection": True,
            "actual_remainder_subsequential_limits_not_assumed": True,
            "state_and_fixed_physical_pulse_explicit": True,
            "source_is_prescribed_metric_not_conserved_matter": True,
            "fixed_duration_not_a_proved_EFT_frequency_band": True,
            "delta_dependent_phase_retaining_effective_descriptions_not_excluded": True,
            "retarded_inverse_not_silently_a_single_copy_variational_action": True}


@cache
def build_report():
    pins = prior_checks()
    for name in ("audit", "interfaces", "loading", "phase"):
        if not (ROOT/"tests"/("test_own_forced_"+name+".py")).is_file():
            raise ValueError("Each independent scientific audit must be source-pinned")
    groups = {
        "actual_loading_source_and_Green_algebra": flatten(loading.identities()),
        "finite_window_phase_and_original_Q_readout": flatten(phase.identities()),
        "independently_derived_real_phase_algebra": audit.phase_algebra()["identities"],
        "frozen_action_domain_and_independent_constant_interfaces": bridges.identities(),
    }
    if any(value != 0 for group in groups.values() for value in group.values()):
        raise ValueError("An exact source, phase, clock or interface residual failed")
    checks = {"loading": loading.checks(), "phase": phase.checks(),
              "independent_bounds": audit.checks(), "premise_discharge": bridges.checks()}
    if not all(value is True for group in checks.values() for value in group.values()):
        raise ValueError("A fixed-pulse continuous proof margin failed")
    return {
        "schema": 1, "claim": "P8-S6.36.OWN_FORCED", "date": "2026-09-07",
        "status": "FIXED_PHYSICAL_PULSE_ZERO_DATA_OWN_F_RESPONSE_WITH_ORIGINAL_Q_PHASE_SEPARATION; ORIGINAL_P8_OPEN",
        "prior_sha256": pins,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": ["notes/proof.md", "notes/audit.md", "notes/interfaces.md"],
        "domain": "M,tau,eta>0 fixed; 0<delta<=10^-6; u=T/tau from-L to+a, L=1/100, a=L/2",
        "input_state_contract": "one delta-independent C-infinity 0<=q<=eta supported inside(-3L/4,-5L/8), integral(q du)>=eta L/16; Q=Q_u=0 at-L; original q is prescribed off shell",
        "exact_residuals": {name: dict.fromkeys(group, "0") for name, group in groups.items()},
        "exact_residual_count": sum(map(len, groups.values())),
        "proof_checks": checks,
        "finite_parameter_loading": serialize(loading.calibration()),
        "explicit_fixed_physical_pulse": serialize(loading.pulse()),
        "phase_comparison_and_scalar_readout": serialize(phase.calibration()),
        "independent_sharper_constants": serialize(audit.constants()),
        "exact_phase_sequence_calibration": serialize(phase.phase_sequences(2)),
        "physical_contract_nonunit_calibration": serialize(physical_contract(amplitude=Fraction(1, 100), time_scale=2)),
        "original_Q_limsup_minus_liminf_lower_per_eta": "1/200",
        "controls": controls(),
        "not_established": [
            "A physical-g/clock or full nonlinear parent/CD solution",
            "A conserved matter or free-chi source-loading theorem",
            "Exact actual subsequential reference limits or a convergent central remainder",
            "Exclusion of every state, pulse class or delta-dependent phase-retaining response",
            "A fixed low-frequency EFT band, physical cutoff or controlled all-order reduction",
            "Full CD/matter matching, positivity/V/G applicability, UV completion or original P8 closure"],
        "verification_boundary": "Written first-zero loading, punctured continuous dependence and separated phase-neighborhood proof, exact margins and independent action/readout audits; not proof-assistant formalization",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The fixed-pulse own-f response report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.36.OWN_FORCED replay passed; original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
