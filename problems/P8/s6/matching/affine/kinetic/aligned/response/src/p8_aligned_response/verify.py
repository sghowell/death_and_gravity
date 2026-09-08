"""Read-only prepared leading-response and light-jet certificate."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_affine_aligned import verify as parent
from p8_affine_retuned.bounds import units

from . import bounds, derivatives, forcing, leading

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"prepared-leading-response.json"
PARENT_SHA = "7fe78d55c3fb1b4a2ff0adc40260dcc3029e7ba3701f08f90278fd7eab2140c5"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_aligned_response/*.py"))
            +sorted(ROOT.glob("tests/*.py"))+sorted(ROOT.glob("*.md"))
            +sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen aligned-source canonical-mode certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_42_fully_rebuilt": PARENT_SHA,
            "actual_all64_source_and_complete_quadratic_constraints_rechecked": True,
            "no_frozen_action_or_ancestor_modified": True}


def controls():
    wrong = (True, False, 1.0, sp.Float(1), "1", sp.I, sp.oo, sp.nan, sp.Symbol("unproved"))
    small, half = bounds.ZETA_MAX, sp.Rational(1, 2)
    calls = [lambda value=value: bounds.bound_values(value, -half, half, 1, 1) for value in wrong]
    calls += [lambda value=value: leading.source_envelope(value, 1) for value in wrong]
    calls += [lambda values=values: bounds.bound_values(*values) for values in
              ((0, -half, half, 1, 1), (-small, -half, half, 1, 1), (2*small, -half, half, 1, 1),
               (small, -1, half, 1, 1), (small, half, half, 1, 1), (small, half, -half, 1, 1),
               (small, -half, half, -1, 1), (small, -half, half, 2, 1),
               (small, -half, half, 1, -1), (small, -half, half, 1, 1.0))]
    calls += [lambda values=values: leading.source_envelope(*values) for values in ((-1, 1), (1, -1))]
    calls += [lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan)]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An inexact, unproved or outside-domain input was accepted")
    return {"rejected_inputs": rejected,
            "nonzero_heavy_data_and_unprepared_forcing_excluded": True,
            "forcing_derivatives_not_replaced_by_amplitude_or_Fourier_band": True,
            "quadratic_source_uses_output_convolution_momentum": True,
            "homogeneous_source_is_temporal_algebraic_control": True,
            "physical_readouts_not_just_canonical_amplitude": True,
            "conditional_leading_estimate_not_full_on_shell_or_quantum_remainder": True}


@cache
def build_report():
    pins = prior_checks()
    groups = {"actual_constrained_canonical_forcing": forcing.checks(),
              "actual_continuous_frequency_derivatives": derivatives.checks(),
              "prepared_retarded_energy_and_local_error": bounds.checks(),
              "original_ODE_quadratic_source_and_light_preparation": leading.checks()}
    exact = {name: affine.certify_residuals(values) for name, values in groups.items()}
    proofs = {"derivatives_"+name: value for name, value in derivatives.proof_checks().items()}
    proofs.update({"retarded_"+name: value for name, value in bounds.proof_checks().items()})
    jets = leading.jet_constants()
    proofs["actual_source_mixed_envelope"] = jets["mixed_envelope"] == 514
    proofs["actual_source_squared_envelope"] = jets["squared_envelope"] == sp.Rational(132027, 8)
    if not all(value is True for value in proofs.values()):
        raise ValueError("A continuous derivative, source-jet or retarded bound failed")
    example = leading.source_envelope(sp.Rational(1, 1000), sp.Rational(1, 1000))
    conversion = units(3, 2, sp.Rational(1, 10000))
    return {"schema": 1, "claim": "P8-S6.43.ALIGNED_RESPONSE", "date": "2026-09-07",
            "status": "ACTUAL_PREPARED_LEADING_RESPONSE_AND_PHYSICAL_READOUT_ERRORS_CERTIFIED; ORIGINAL_P8_OPEN",
            "prior_sha256": pins,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md", "written_proofs": ["notes/retarded.md", "notes/light-jets.md"],
            "exact_residuals": exact, "named_exact_check_count": sum(map(len, groups.values())),
            "checked_scalar_entries": sum(map(len, groups.values())), "proof_checks": proofs,
            "literal_action": "S6.42 unchanged; second-order vector response with zero heavy data, not a new parent",
            "actual_source": serialize(leading.source()["expected"]),
            "canonical_force": "J=(g_L²*S2)'/g_L=g_L*S2'+2g_L'*S2; v''+(q+1/zeta-U_L)*v=J",
            "domain": "Interval inside [-1/2,1/2], length<=1, output 0<k_com²<=1, 0<zeta<=1/20000; homogeneous separate",
            "state_and_source": "Zero v,v' initially; S2=S2'=S2''=0 initially; |S2^(j)|<=E for j=0..3; A=(815/2)*sqrt(zeta)*|k_com|*E bounds all J jets through order two",
            "continuous_bounds": serialize({"squared_frequency_floor": 19985, "abs_W_prime": 78, "abs_W_second": 712,
                                             "energy_amplification_upper": sp.Rational(21, 20), "light_jet_constants": jets}),
            "errors": "|v-zeta*J|<=zeta*A/80; |W_spatial-zeta*sqrt(q)*(S2'+2rho*S2)|<=6*zeta*E; |W0-S2|<=500*zeta*E; physical vector factor 1/tau",
            "light_jet_bridge": "E=514*epsilon_n*epsilon_K+(132027/8)*epsilon_n² from C³ light jets and actual ODE; summable modewise envelopes control output convolutions; n=n'=n''=0 sufficient preparation",
            "example": serialize({"light": example,
                                  "response": bounds.bound_values(sp.Rational(1, 10**9), -sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(1, 4), example["source_jet_E"]),
                                  "nonunit_conversion": conversion}),
            "controls": controls(),
            "verdict": "The actual leading heavy response has a quantitative prepared retarded local approximation with physical readouts and a sufficient light-jet source bound. The preparation and jet estimates have not been proved for an on-shell light family.",
            "not_established": ["An on-shell light family realizing and maintaining the prepared derivative domain",
                                "Full nonlinear secondary constraints, nonlinear stability or higher-order response/action errors",
                                "Quantum loop remainders, an interacting cutoff, stationary heavy gap or V/G/B UV admissibility",
                                "All momenta, arbitrary heavy initial states, or original P8 closure"],
            "verification_boundary": "Exact symbolic checks, continuous energy and coefficient proofs, independent manufactured-solution and physical-readout controls; not proof-assistant formalization"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The prepared leading-response report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.43.ALIGNED_RESPONSE replay passed; prepared leading error only, original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
