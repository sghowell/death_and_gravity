"""Read-only actual finite-width reference quadratic bound, with exact errors."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as certificate
from p8_proca_relational_qei import verify as parent

from . import audit, bounds, canonical, compact, intervals, proofs, time_jets

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"actual-finite-width-relational-bound.json"
PARENT_SHA="e9ae3ec1449071bddb04cb3106fc1f82194d8b885f6357a39405849d5a5df5d0"
sha,serialize=certificate.sha,certificate.serialize

def source_files():
    return (sorted(ROOT.glob("src/p8_proca_finite_width/*.py"))+sorted(ROOT.glob("tests/*.py"))
        +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))

@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen relational quadratic-inequality report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_100_fully_rebuilt":PARENT_SHA,
        "same_original_global_action_density_source_observable_and_state_class":True,
        "new_quantitatively_prepared_reference_state_not_a_free_matter_replacement":True}

@cache
def build_report():
    prior=prior_checks()
    identities,gates=audit.residuals(),audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("An actual finite-width inequality proof gate failed")
    exact=certificate.certify_residuals(identities)
    jet=intervals.build()
    jetdata={key:value for key,value in jet.items()
        if key.endswith("integers") or key.startswith("frequency_lower")}
    jetdata.update({"working_precision_bits":160,"whole_real_interval":["-1/100","1/100"],
        "Riccati_order":6,"time_jet_order":9,"independent_exact_center_inclusions":proofs.exact_jet_inclusions()})
    initial={"bounds":compact.center_bounds(),"center_Riccati_coefficients":time_jets.build(8)["at_center"],
        "Borel_initial_tail_relative_to_sixth_order":"norm <= 600001/k^8 for k>=16",
        "infrared_covariance_mixing_interval":[16,32],"regular_initial_covariance_norm_without_hbar_over_two_kappa":16,
        "fresh_global_generalized_Hadamard_state_family":True}
    return {"schema":1,"claim":"P8-S6.101.ACTUAL_FINITE_WIDTH_RELATIONAL_BOUND","date":"2026-09-09",
        "status":"COMPLETE_EXPLICIT_NONOPTIMAL_ACTUAL_GLOBAL_RELATIONAL_TEST_ENERGY_REFERENCE_DIFFERENCE_BOUND_FOR_ALL_COMPACT_COMOVING_SAMPLERS_IN_ABS_U_LT_ONE_HUNDREDTH_WITH_POSITIVE_ALL_ORDER_INITIAL_STATE_WHOLE_INTERVAL_REMAINDERS_ZERO_MOMENTUM_AND_COMPLETE_MODE_ERROR; PHYSICAL_STRESS_SHARP_COSMOLOGICAL_CALIBRATION_V_G_B_AND_ORIGINAL_P8_OPEN",
        "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
        "formulation":"FORMULATION.md","written_proofs":["notes/model.md","notes/state.md","notes/remainder.md",
            "notes/fourier.md","notes/bound.md","notes/scope.md"],
        "exact_residuals":exact,"named_exact_check_count":len(identities),
        "checked_scalar_entries":sum(value.rows*value.cols if isinstance(value,sp.MatrixBase) else 1
            for value in identities.values()),"proof_checks":gates,
        "actual_all_real_frequency_chart":serialize(canonical.data()),
        "actual_initial_state_and_eighth_order_coefficients":serialize(initial),
        "whole_interval_jet_and_remainder_bounds":serialize(jetdata),
        "explicit_actual_finite_width_inequality":serialize(bounds.data()),
        "controls":audit.controls(),
        "verdict":"For a nonempty explicitly constrained family of actual globally transported generalized Hadamard reference states and every real smooth comoving sampler supported in |u|<1/100, the original relational test-energy reference functional is bounded by fully specified positive constants multiplying ||g||_2^2 and two squared linear combinations of ||g^(j)||_1 for j<=3. The resulting reference-difference inequality retains zero spatial momentum, an all-order positive initial state, the full time-dependent sixth-order graph defect and its nonzero initial tail. The bound is nonoptimal, not a sharp cosmological calibration or the complete covariant physical stress. Original P8 remains open.",
        "not_established":["An optimal bound or controlled relative correction to the leading short-sampling coefficient",
            "Sharp cosmological calibration, a physical-stress QEI or a realistic-field singularity theorem",
            "Covariant renormalized complete physical stress, metric/clock response, Ward identities and self-consistent semiclassical backreaction",
            "Full quantum tensor/vector construction, interacting or omitted-operator/loop/heavy-threshold validity, cutoff or actual vacuum/gravity/common-parent V/G/B matching",
            "Original realistic/null/nonminimal P8(a), original UV classification P8(b) or original P8 closure"],
        "verification_boundary":"Native complete regular canonical and scaled symplectic identities, exact zero-momentum/source/observable bridges, actual Hamiltonian parity, rational compact energy and root-isolated eighth-order state bounds, all triangular exact Riccati residuals through eighth and ninth order, 196 independent exact jet inclusions, whole-interval 160-bit ball automatic differentiation, explicit positive-frequency resolvent and ordered IBP coefficient recurrences, and exact low/high/error integrals and normalization. Written all-order Borel preparation, positive global state, Duhamel mode remainder, matrix Fourier estimate and reference-difference proofs are pinned, not proof-assistant formalized. Every ancestor, source/proof/test byte and report field is independently replayed. Scientific SymPy is unchanged."}

def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The actual finite-width relational-bound report differs from read-only replay")

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.101.ACTUAL_FINITE_WIDTH_RELATIONAL_BOUND replay passed; full physical stress, V/G/B and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))

if __name__=="__main__":
    main()
