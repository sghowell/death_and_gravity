"""Read-only exact nearby classical retarded relational source certificate."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as original
from p8_proca_nearby_packets import verify as parent

from . import audit, fronts, remainder, source, support, transport

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"retarded-nearby-relational-response.json"
PARENT_SHA="5f840f8eae3eeea353236e4106743bc412529d61948c9f0522dfba9564dfa01b"
sha,serialize=original.sha,original.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_proca_nearby_retarded/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen finite-frequency nearby classical packet report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_90_fully_rebuilt":PARENT_SHA,
            "same_actual_nearby_classical_solution_and_literal_source_response":True,
            "no_quantum_state_profile_or_interacting_EFT_cutoff_transfer":True}


@cache
def build_report():
    prior=prior_checks()
    identities,gates=audit.residuals(),audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("An actual classical retarded-source proof gate failed")
    exact=original.certify_residuals(identities)
    return {"schema":1,"claim":"P8-S6.91.CLASSICAL_NEARBY_RETARDED","date":"2026-09-09",
            "status":"EXACT_NEARBY_CLASSICAL_LOCAL_RELATIONAL_RETARDED_RESPONSE_HAS_A_NONZERO_CLOCK_FRONT_OUTSIDE_THE_PHYSICAL_MATTER_CONE_AND_NONZERO_COMPACT_PROBE_PAIRINGS; FINITE_EFT_QUANTUM_MATCHING_AND_ORIGINAL_P8_OPEN",
            "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
            "formulation":"FORMULATION.md",
            "written_proofs":["notes/source.md","notes/transfer.md","notes/distributions.md","notes/compact.md","notes/scope.md"],
            "exact_residuals":exact,"named_exact_check_count":len(identities),
            "checked_scalar_entries":sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in identities.values()),
            "proof_checks":gates,
            "literal_local_relational_source_and_complete_canonical_maps":serialize(source.data()),
            "uniform_two_time_transfer_with_both_endpoint_changes":
                serialize({"constants":transport.constants(),
                           "sample_errors":{str(k):transport.error(k) for k in
                                            (transport.normal_form.MIN_MOMENTUM,sp.Integer(10)**72)}}),
            "exact_clock_and_matter_front_amplitudes":serialize(fronts.data()),
            "complete_frequency_Fourier_and_L2_remainder_bounds":serialize(remainder.data()),
            "actual_matter_spacelike_compact_probe_neighborhoods":serialize(support.data()),
            "controls":audit.controls(),
            "verdict":"The exact linear classical response to sqrt(-g) J[chi-chi_background(phi)] on the actual nearby S6.88 bounce is the sum of positive clock and matter delta shells and an L2 spatial remainder, uniformly for ordered inner-interval times. The full original Hamiltonian controls all compact low frequencies, and both endpoint normal-form changes give a uniform high-frequency multiplier error below 10^44/k^2. The nonzero clock shell cannot be cancelled by the L2 remainder or the separate matter shell. Some real smooth compact source and detector probes, with every support pair matter-spacelike, have nonzero retarded response. The front uses unbounded classical frequencies and is not transferred to a finite interacting EFT, UV parent or quantum commutator.",
            "not_established":["A specified compact pulse or detector with quantitative nonzero signal, noise, resolution or finite admitted frequency band",
                "A classical-to-interacting-EFT front/response transfer, canonical interaction cutoff or omitted-operator, loop and heavy-threshold error",
                "A quantum state, fixed profile, cancellation, response or inverse on these nearby backgrounds",
                "A nonlinear source apparatus, nonlinear sourced solution, all-wavelength numerical stability or globally complete nearby bounce",
                "Finite Wilson matching, healthy full-candidate vacuum, common-parent V/G/B, UV completion/exclusion or original P8 closure"],
            "verification_boundary":"Exact native original-chart source and moving canonical identities, both-time endpoint transfer algebra, actual mode-pair amplitudes and bounds, full compact-frequency matrix majorants, Fourier and Plancherel integrals, exact open-neighborhood margins and written distribution/product-test proofs. Every parent is natively rebuilt; no floating inequality, scientific-library patch, discarded frequency range, quantum transfer or cutoff assertion. Not proof-assistant formalization."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The nearby classical retarded relational response report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.91.CLASSICAL_NEARBY_RETARDED replay passed; finite EFT, quantum, matching and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
