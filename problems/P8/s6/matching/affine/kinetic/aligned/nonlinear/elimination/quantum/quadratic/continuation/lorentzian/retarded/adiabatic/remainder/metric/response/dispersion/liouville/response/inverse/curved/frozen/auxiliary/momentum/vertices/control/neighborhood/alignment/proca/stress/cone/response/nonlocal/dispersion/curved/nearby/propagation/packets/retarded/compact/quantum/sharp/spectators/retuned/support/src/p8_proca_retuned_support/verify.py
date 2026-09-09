"""Read-only complete retuned scalar support and new quadratic CCR certificate."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as certificate
from p8_proca_retuned_margin import verify as parent

from . import audit, bridge, domain, growth, momentum, quantum

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"retuned-scalar-support-microcausality.json"
PARENT_SHA="6e41de7d01b5685b98a6f6659ff9c7455af9fa98b0e60e0b96b7ef70cfd61e03"
sha,serialize=certificate.sha,certificate.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_proca_retuned_support/*.py"))+sorted(ROOT.glob("tests/*.py"))
        +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen new covariant-margin report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_96_fully_rebuilt":PARENT_SHA,
        "support_and_state_use_new_literal_action_and_new_actual_solution":True,
        "old_fast_clock_obstruction_and_old_report_bytes_remain_distinct":True}


@cache
def build_report():
    prior=prior_checks()
    identities,gates=audit.residuals(),audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A complete new-support or quadratic-state gate failed")
    exact=certificate.certify_residuals(identities)
    return {"schema":1,"claim":"P8-S6.97.RETUNED_SCALAR_SUPPORT_MICROCAUSALITY","date":"2026-09-09",
        "status":"COMPLETE_RETUNED_LOCAL_SCALAR_RESPONSE_HAS_ACTUAL_MATTER_CAUSAL_SUPPORT_AND_FRESH_REDUCED_QUADRATIC_CCR_STATE_IS_RELATIONALLY_MICROCAUSAL_ON_CERTIFIED_INTERVAL; INTERACTING_QUANTUM_V_G_B_AND_ORIGINAL_P8_OPEN",
        "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
        "formulation":"FORMULATION.md",
        "written_proofs":["notes/domains.md","notes/evolution.md","notes/support.md","notes/quantum.md","notes/scope.md"],
        "exact_residuals":exact,"named_exact_check_count":len(identities),
        "checked_scalar_entries":sum(value.rows*value.cols if isinstance(value,sp.MatrixBase) else 1
                                     for value in identities.values()),
        "proof_checks":gates,
        "new_literal_complete_source_and_frequency_bridge":serialize(bridge.data()),
        "fresh_actual_coefficient_and_holomorphic_basis_bounds":serialize({
            "actual":domain.actual_bounds(),"basis":domain.basis(),"complete":domain.coefficient_majorants()}),
        "complex_spatial_momentum_geometry_and_null_cone_coverage":serialize(momentum.data()),
        "complete_complex_frequency_growth_and_actual_matter_radius":serialize(growth.data()),
        "new_quadratic_density_state_and_causal_Kubo_response":serialize(quantum.data()),
        "controls":audit.controls(),
        "verdict":"For the separately named epsilon=1/200 action and its own certified nearby solution, the complete original local relational scalar response is entire in complex spatial momentum. Every inverse-frequency term and the original complex-null-cone chart are retained. Fresh logarithmic-norm estimates give global exponential type equal to the actual physical-matter null distance, and the written distributional support argument gives the complete retarded matter-cone support on [-5*10^-8,5*10^-8]. A fresh positive tempered Gaussian state transported with the new density-canonical evolution has the exact new Kubo commutator, hence quadratic relational microcausality on this interval. This is not an interacting, Hadamard, semiclassical or UV-completion result.",
        "not_established":["A common quantitative full-support interval for the entire central family or global nearby continuation",
            "Nonlinear stability, canonical interactions, omitted operators, loops or an interacting EFT validity band",
            "Hadamard state, renormalized stress, self-consistent semiclassical backreaction or interacting microcausality",
            "Physical vacuum, amplitudes, dispersion and common-parent V/G/B matching",
            "Original realistic-field P8(a), original UV classification P8(b) or original P8 closure"],
        "verification_boundary":"Native complete source/Laurent and density-canonical identities, fresh rational whole-box coefficient bounds and holomorphic basis majorants, exact logarithmic-norm principal cancellation, complex root Gram identity and explicit null-vector coverage, positive covariance factorization and Kubo normalization. Pinned written Peano-Baker entire-dependence, complete matrix Gronwall, exact-radius distributional support, Schwartz/Duhamel and Gaussian-state arguments. All ancestors are natively rebuilt and every source/proof/test byte pinned; every report field mutation is rejected. No numerical point-scan proof, scientific-library patch, old state transfer or proof-assistant formalization."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The retuned scalar-support report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.97.RETUNED_SCALAR_SUPPORT_MICROCAUSALITY replay passed; interacting quantum V/G/B and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
