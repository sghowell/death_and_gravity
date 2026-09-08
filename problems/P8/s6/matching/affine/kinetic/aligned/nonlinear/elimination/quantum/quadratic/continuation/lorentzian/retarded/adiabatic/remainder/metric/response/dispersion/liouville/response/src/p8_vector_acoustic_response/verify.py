"""Read-only prepared acoustic covariance and physical stress-response bridge."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_vector_liouville import verify as parent

from . import dimensional, proofs, readout, source, transport

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"acoustic-covariance-response.json"
PARENT_SHA = "daab6d10f0dfa4190153fa63f5bb6eed9829c1ae73a0704eac34efb5c5b0da10"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_vector_acoustic_response/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen acoustic potential-reduction report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_70_fully_rebuilt": PARENT_SHA, "physical_state_and_finite_prescription_unchanged": True}


def controls():
    calls = []
    for value in (True, False, sp.true, sp.false, 0, 1, 1.0, sp.Integer(1), "transverse", "", None):
        calls.extend((lambda value=value: source.data(value),
                      lambda value=value: source.derivative(source.variation.n[0], value)))
    for value in (True, False, sp.true, sp.false, 1.0, sp.Float(1), "0", None,
                  source.variation.n[0]*source.variation.zeta[0], source.variation.n[0]*source.xi_L,
                  sp.Integer(1), sp.Symbol("foreign")*source.variation.n[0], sp.oo, sp.nan):
        calls.append(lambda value=value: source.clean(value))
    calls.extend(lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan))
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An unsupported acoustic response input was accepted")
    return {"rejected_inputs": rejected, "native_validation_before_cached_implementations": True,
            "independent_prepared_clock_shifts_retained": True,
            "nonzero_initial_squeezing_variance_kernel_fixture": "-3/2",
            "state_independent_Green_does_not_make_variance_response_state_independent": True,
            "physical_output_and_second_mass_contacts_retained": True,
            "no_early_dimensional_transverse_pump_deletion": True,
            "no_new_auxiliary_scalar_subtraction_or_C2_stress_bound": True}


@cache
def build_report():
    prior = prior_checks()
    identities, gates = proofs.residuals(), proofs.checks()
    if not all(value is True for value in gates.values()):
        raise ValueError("An acoustic covariance-response audit gate failed")
    certified = affine.certify_residuals(identities)
    return {"schema": 1, "claim": "P8-S6.71.ACOUSTIC_COVARIANCE_RESPONSE", "date": "2026-09-08",
            "status": "EXACT_PREPARED_ACOUSTIC_COVARIANCE_AND_PHYSICAL_READOUT_BRIDGE; ORIGINAL_P8_OPEN",
            "prior_sha256": prior, "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
            "formulation": "FORMULATION.md", "written_proofs": ["notes/transport.md", "notes/readout.md", "notes/limit.md"],
            "exact_residuals": certified, "named_exact_check_count": len(identities),
            "checked_scalar_entries": sum(v.rows*v.cols if isinstance(v, sp.MatrixBase) else 1 for v in identities.values()),
            "proof_checks": gates,
            "generic_covariance_transport": serialize({"physical_to_acoustic": transport.covariance_map(),
                "acoustic_generator": transport.matrix(transport.bk, transport.D),
                "varied_acoustic_generator": transport.varied_generator()}),
            "prepared_sector_covariance_sources": serialize({sector: source.data(sector) for sector in ("T", "L")}),
            "physical_readout_response_and_contact_rows": serialize(readout.actual_readouts()),
            "generic_physical_readout_variation": serialize({key: value for key, value in readout.algebra().items() if key != "checks"}),
            "full_dimensional_pump_squares": serialize(dimensional.pumps()),
            "retarded_initial_policy": "The exact selected normalized modes define the covariance response. Varied initial data vanish on the prepared neighborhood; unvaried initial squeezing is retained. The two acoustic shifts start at zero but need not vanish at the final endpoint.",
            "finite_limit_policy": "Use the same physical adiabatic subtraction mapped through the positive D-independent clocks, retaining D-1 transverse multiplicity, both dimensional pump jets, analytic plus/minus modes and the S6.67 finite local term. The S6.68 finite C10-to-C0 response is unchanged; no independent scalar subtraction is introduced.",
            "controls": controls(),
            "verdict": "Exact moving-clock covariance and tangent transport reproduce both physical source vertices and all four physical energy/pressure response rows. The retarded mode Green function has the fixed normalization and Kubo sign, while its variance kernel retains the selected state. Both prepared time shifts, physical output/mass contacts and dimensional pump jets remain. This is an exact representation of the already finite prepared response, not a new no-loss stress or coupled-inverse theorem.",
            "not_established": ["A C2-to-C0 norm for the full stress, an integrated logarithmic-remainder bound or a coupled causal inverse",
                                "Spatial response, independent initial states, quantum stability/cones, interactions, other loops, cutoff, finite Wilson matching, V/G/B or original P8 closure"],
            "verification_boundary": "Exact covariance, Green, source, physical readout and dimensional identities with written prepared-state and unchanged-regulator arguments. Not proof-assistant formalization."}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The acoustic covariance-response report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.71.ACOUSTIC_COVARIANCE_RESPONSE replay passed; coupled inverse and original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
