"""Read-only exact finite-frequency classical nearby relational packets."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as original
from p8_proca_nearby_cones import verify as parent

from . import audit, domains, finite, modes, normal_form, observable, packets

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"finite-frequency-nearby-packets.json"
PARENT_SHA="f176876ff1cbc122fc9cc37c4b63cf8b1494f81ec73ddad53f2771f56e7af530"
sha,serialize=original.sha,original.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_proca_nearby_packets/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen actual nearby classical cone report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_89_fully_rebuilt":PARENT_SHA,
            "same_actual_nearby_classical_solution_and_literal_linearized_action":True,
            "no_quantum_state_profile_or_interacting_band_transfer":True}


@cache
def build_report():
    prior=prior_checks()
    identities,gates=audit.residuals(),audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("An actual finite-frequency nearby packet gate failed")
    exact=original.certify_residuals(identities)
    return {"schema":1,"claim":"P8-S6.90.CLASSICAL_NEARBY_PACKETS","date":"2026-09-09",
            "status":"EXACT_FINITE_FREQUENCY_CLASSICAL_NEARBY_RELATIONAL_MATTER_PACKETS_HAVE_CONTROLLED_MODE_DIFFRACTION_AND_TAIL_ERRORS_AND_ADVANCE_RELATIVE_TO_MATTER_CONE; COMPACT_PREPARATION_INTERACTING_CUTOFF_QUANTUM_MATCHING_AND_ORIGINAL_P8_OPEN",
            "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
            "formulation":"FORMULATION.md",
            "written_proofs":["notes/canonical.md","notes/analytic.md","notes/observable.md","notes/packets.md","notes/scope.md"],
            "exact_residuals":exact,"named_exact_check_count":len(identities),
            "checked_scalar_entries":sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in identities.values()),
            "proof_checks":gates,
            "exact_finite_frequency_canonical_system":serialize(finite.data()),
            "action_normalized_modes_and_complete_transport":serialize(modes.data()),
            "actual_complex_solution_domain_and_annuli":serialize(domains.analytic_domain()),
            "literal_finite_frequency_Laurent_coefficient_majorants":serialize(domains.coefficient_majorants()),
            "explicit_analytic_basis_and_inverse_majorants":serialize(domains.basis_majorants()),
            "actual_near_identity_normalform_and_uniform_constants":
                serialize({"matrices":normal_form.data(),"constants":normal_form.constants()}),
            "local_relational_observable_pole_and_reconstruction":serialize(observable.data()),
            "exact_three_dimensional_classical_packets_and_normalized_errors":
                serialize({str(value):packets.packet(value) for value in (72,75,120)}),
            "controls":audit.controls(),
            "verdict":"On the inner real half of the actual S6.88/S6.89 nearby classical bounce interval, exact scalar linearized solutions form outgoing three-dimensional real packets in O=chi-chi_background(phi). The full finite-frequency Hamiltonian, volume rescaling and symplectic mode basis retain every contact and inverse-frequency remainder. Leading diagonal transport vanishes exactly. Actual complex-time bounds give a quantitative normal form. At comoving carrier 10^72 and bandwidth 10^24, modal reconstruction plus diffraction error is below 10^-27 and normalized exact outside-cube L2 mass is below 10^-8. The clock comparison cube lies outside the matter causal future of the initial comparison cube, with nonzero packet tails retained. The result is finite-frequency classical transport, not compact local preparation, a quantum signal or a frequency band admitted below an interacting EFT cutoff.",
            "not_established":["Compact spatial support, compact local retarded preparation or a quantum detector/commutator theorem",
                "A classical-to-interacting-EFT frequency transfer, optimized carrier, canonical interaction control, omitted operators, loops or cutoff",
                "A quantum state, fixed profile, cancellation, response or inverse on these nearby backgrounds",
                "Exact nonlinear packet solutions, all-wavelength stability or globally complete nearby bounces",
                "Heavy thresholds, finite Wilson matching, healthy full-candidate vacuum, common-parent V/G/B, UV completion/exclusion or original P8 closure"],
            "verification_boundary":"Exact full finite-frequency Hamiltonian identities, symplectic/action normalization and zero diagonal transport, actual complex solution annuli and complete Laurent bounds, time Cauchy derivatives, exact near-identity remainder and rational error constants, relational gauge/pole checks, explicit Fourier normalization, diffraction and normalized nonzero-tail estimates, with written ODE and localization proofs. Every parent is natively rebuilt. No floating inequality, library patch, quantum transfer, compact-support assumption or cutoff assertion. Not proof-assistant formalization."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The finite-frequency nearby packet report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.90.CLASSICAL_NEARBY_PACKETS replay passed; compact preparation, EFT cutoff, quantum, matching and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
