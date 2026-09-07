"""Read-only exact own-f fixed-physical-window transfer certificate."""
import argparse
import hashlib
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

import sympy as sp
from p8_exact_stationary import intervals
from p8_exact_stationary import verify as prior

from . import action, arb_audit, bridges, connection, connection_audit, potential

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"own-fixed-window-transfer.json"
PARENT_SHA = "727811a8563584cc511e749ff038eadaa416f02aae268a15fa4f90eb8ca4b6b3"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_files():
    return (sorted(ROOT.glob("src/p8_own_scattering/*.py"))
            +sorted(ROOT.glob("tests/*.py"))+sorted(ROOT.glob("*.md"))
            +sorted(ROOT.glob("notes/*.md")))


def serialize(value):
    """Exact real/complex formal algebra, rational intervals, or JSON metadata.

    Unlike a physical real-valued input, a named wave-basis coefficient is
    necessarily complex. Such symbolic output is allowed explicitly; no
    binary float or nonfinite expression is admitted as certified data.
    """
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, intervals.Interval):
        return {"lower": str(value.lo), "upper": str(value.hi)}
    if isinstance(value, sp.MatrixBase):
        return serialize(value.tolist())
    if isinstance(value, sp.Basic):
        if value.has(sp.Float, sp.oo, -sp.oo, sp.zoo, sp.nan):
            raise TypeError("Only exact finite formal report algebra is admitted")
        return str(value)
    if isinstance(value, dict):
        return {str(key): serialize(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(item) for item in value]
    if value is None or isinstance(value, (str, int, bool)):
        return value
    raise TypeError("Only exact algebra or JSON metadata is admitted")


@cache
def prior_checks():
    if sha(prior.REPORT) != PARENT_SHA:
        raise ValueError("The frozen exact own-f stationary branch changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {"S6_33_recursively_rebuilt": PARENT_SHA,
            "same_unchanged_action_metric_clock_and_free_chi": True,
            "sibling_of_S6_34_not_dependent_on_its_shrinking_pulse": True}


def physical_contract(delta=Fraction(1, 10**6), time_scale=1):
    delta, tau = map(intervals.rational, (delta, time_scale))
    if not 0 < delta <= connection.DELTA_MAX or tau <= 0:
        raise ValueError("Require 0<delta<=10^-6 and positive physical time scale")
    return {"delta": delta, "tau": tau,
            "physical_half_width": tau*connection.HALF_WIDTH,
            "physical_duration": 2*tau*connection.HALF_WIDTH,
            "physical_potential_remainder_upper": potential.POTENTIAL_CAP/tau**2,
            "dimensionless_velocity": "Q_u=tau*Q_T",
            "endpoint_wave_state": "(psi,psi_t/rho), with actual k and k_u",
            "homogeneous_mixing_lower": Fraction(1, 80),
            "fixed_temporal_band_or_EFT_cutoff": False,
            "zero_data_matter_response": False}


def controls():
    bad = (True, False, 1.0, sp.Float("1"), sp.oo, sp.nan, "1/1000000", sp.Symbol("delta"))
    calls = [lambda value=value: physical_contract(delta=value) for value in bad]
    calls += [lambda value=value: physical_contract(delta=value)
              for value in (0, -1, Fraction(1, 100))]
    calls += [lambda value=value: physical_contract(time_scale=value)
              for value in (True, 1.0, 0, -1)]
    calls += [lambda value=value: serialize(value)
              for value in (0.1, sp.Float("0.1"), sp.oo, -sp.oo, sp.zoo, sp.nan)]
    calls += [lambda value=value: connection.endpoint_map(0, Fraction(1, 10**6), value, 0)
              for value in (True, sp.true, 0.1, 0, -1, sp.I)]
    calls += [lambda: potential.center_data(0),
              lambda: potential.center_data(0, extension=1),
              lambda: potential.partial_jets(root_bracket=intervals.Interval(3)),
              lambda: arb_audit.enclosure(64),
              lambda: arb_audit.enclosure(True)]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An invalid arithmetic, root, domain or physical input was accepted")
    return {"rejected_inputs": rejected,
            "literal_delta_zero_not_an_action": True,
            "bounded_remainder_not_jointly_continuous_or_derivative_bounded": True,
            "own_f16_not_coupled_relative80": True,
            "actual_k_u_and_physical_time_column_retained": True,
            "ordinary_Gauss_not_Olver_normalized_function": True,
            "time_transfer_det_plus_one_radial_det_minus_one": True,
            "full_reference_tails_included": True,
            "conditional_remainder_premise_discharged_in_this_report": True,
            "fictitious_free_exterior_not_physical_vacuum": True,
            "fixed_physical_window_not_fixed_low_frequency_EFT_band": True,
            "homogeneous_data_not_zero_data_source_loading": True,
            "retarded_single_copy_variational_action_not_assumed": True}


@cache
def build_report():
    pins = prior_checks()
    for name in ("arb_audit", "connection", "connection_audit", "interfaces", "jets", "potential"):
        if not (ROOT/"tests"/("test_own_scattering_"+name+".py")).is_file():
            raise ValueError("Every separately authored scientific audit must be pinned")
    groups = {"independent_canonical_action_source_and_clock": action.identities(),
              "reference_connection_and_finite_window_algebra": connection.identities(),
              "frozen_action_profile_endpoint_and_domain_interfaces": bridges.identities()}
    if any(value != 0 for group in groups.values() for value in group.values()):
        raise ValueError("An exact action, connection or interface residual failed")
    proof_checks = {"whole_box_potential": potential.checks(),
                    "reference_and_transfer_margins": connection.checks(),
                    "premise_discharge": bridges.checks()}
    if not all(value is True for group in proof_checks.values() for value in group.values()):
        raise ValueError("A continuous-domain bound or conditional interface failed")
    connection_balls = connection_audit.report()
    if not all(value is True for value in connection_balls["checks"].values()):
        raise ValueError("An independent Acb normalization or strict-margin audit failed")
    return {
        "schema": 1, "claim": "P8-S6.35.OWN_SCATTERING", "date": "2026-09-07",
        "status": "EXACT_OWN_F_FIXED_PHYSICAL_WINDOW_HOMOGENEOUS_TRANSFER; ORIGINAL_P8_OPEN",
        "prior_sha256": pins,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": ["notes/potential.md", "notes/connection.md", "notes/interfaces.md"],
        "source_audit": "notes/sources.md",
        "potential_domain": "|u|<=1/100, 0<delta<=1/100; u=T/tau, c=2+delta",
        "transfer_domain": "|T|<=tau/100, 0<delta<=10^-6; M,tau>0",
        "probe_and_state_contract": "same off-shell physical g=eta and clock theta=T; q=0, prescribed incoming homogeneous Q data; exact own-f isotropic background only",
        "exact_residuals": {name: dict.fromkeys(group, "0") for name, group in groups.items()},
        "exact_residual_count": sum(map(len, groups.values())),
        "proof_checks": proof_checks,
        "fraction_partial_jet_potential_proof": serialize(potential.calibration()),
        "independent_Arb_total_jet_cover": arb_audit.report(),
        "independent_Acb_connection_corroboration": connection_balls,
        "canonical_action_source_and_endpoint_maps": serialize(action.calibration()),
        "reference_Gauss_connection": serialize(connection.coefficients()),
        "reference_parity_basis": serialize(connection.even_odd()),
        "reference_time_transfer": serialize(connection.time_transfer()),
        "finite_window_margins": serialize(connection.calibration()),
        "physical_contract_nonunit_calibration": serialize(physical_contract(time_scale=2)),
        "controls": controls(),
        "not_established": [
            "A physical-g/clock or full nonlinear parent/CD solution",
            "A fixed-source zero-data response or conserved matter/free-chi source matching",
            "An actual asymptotic parent vacuum, quantum state or physical scattering S matrix",
            "A vanishing-error outer dressing, uniform derivative expansion or all-order remainder",
            "A fixed physical low-frequency EFT band, cutoff, scalar health or full CD dictionary",
            "Positivity/V/G applicability, UV completion or original P8 closure"],
        "verification_boundary": "Written exact-jet, pole-cancellation, source-audited connection and energy/Duhamel proof; independent Fraction and Arb derivations, Acb corroboration; not proof-assistant formalization",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The own-f fixed-window report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.35.OWN_SCATTERING replay passed; original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
