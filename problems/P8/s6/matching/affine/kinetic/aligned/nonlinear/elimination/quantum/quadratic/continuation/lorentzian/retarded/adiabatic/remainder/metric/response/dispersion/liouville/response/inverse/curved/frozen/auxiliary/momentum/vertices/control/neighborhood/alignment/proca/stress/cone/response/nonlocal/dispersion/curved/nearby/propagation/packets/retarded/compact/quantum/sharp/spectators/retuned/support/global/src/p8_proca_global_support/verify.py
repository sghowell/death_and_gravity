"""Read-only complete global retuned scalar support and new quadratic CCR certificate."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as certificate
from p8_proca_retuned_support import verify as parent

from . import audit, charts, coverage, model, modes, phase, quantum, transfer

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"global-retuned-scalar-support.json"
PARENT_SHA="128686063857028a8e8d1b1f00818e311811c43440d4bcbfde0b033fabaf2b35"
sha,serialize=certificate.sha,certificate.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_proca_global_support/*.py"))+sorted(ROOT.glob("tests/*.py"))
        +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen nearby retuned scalar-support report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_97_fully_rebuilt":PARENT_SHA,
        "same_named_retuned_action_but_separate_global_clock_solution":True,
        "global_state_not_transferred_from_either_nearby_solution":True}


@cache
def build_report():
    prior=prior_checks()
    identities,gates=audit.residuals(),audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A complete global-support chart or pole-separation gate failed")
    exact=certificate.certify_residuals(identities)
    return {"schema":1,"claim":"P8-S6.98.GLOBAL_RETUNED_SCALAR_SUPPORT","date":"2026-09-09",
        "status":"COMPLETE_SCALAR_RESPONSE_ON_PRESERVED_GLOBAL_RETUNED_CLOCK_HAS_ACTUAL_MATTER_CAUSAL_SUPPORT_AT_ALL_FINITE_TIMES_AND_NEW_GLOBAL_REDUCED_QUADRATIC_CCR_STATE_IS_RELATIONALLY_MICROCAUSAL; INTERACTING_QUANTUM_V_G_B_AND_ORIGINAL_P8_OPEN",
        "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
        "formulation":"FORMULATION.md",
        "written_proofs":["notes/canonical.md","notes/charts.md","notes/growth.md","notes/quantum.md","notes/scope.md"],
        "exact_residuals":exact,"named_exact_check_count":len(identities),
        "checked_scalar_entries":sum(value.rows*value.cols if isinstance(value,sp.MatrixBase) else 1
                                     for value in identities.values()),
        "proof_checks":gates,
        "actual_retuned_global_clock_background_and_coefficient_substitution":serialize({
            "background":model.background(),"actual_coefficient_substitution":model.substitution()}),
        "complete_regular_density_phase_and_original_local_source":serialize(phase.data()),
        "complete_unitary_and_gamma_Euler_systems":serialize(charts.data()),
        "full_chart_intertwining_and_polynomial_endpoint_weights":serialize(transfer.data()),
        "explicit_global_chart_cover_and_complex_pole_separation":serialize(coverage.data()),
        "exact_principal_modes_and_real_frequency_growth_cancellation":serialize(modes.data()),
        "new_global_quadratic_scalar_state_and_causal_Kubo_response":serialize(quantum.data()),
        "controls":audit.controls(),
        "verdict":"For the separately named total-margin epsilon=1/200 action, the preserved original global classical clock a=(1+u^2)^2 has complete reduced local scalar retarded support in its actual physical-matter cone for every finite source/output time. Both complete Euler charts are exactly intertwined with the regular density phase; all finite-q velocity poles and complex spatial null directions are covered. Explicit chart overlaps and complex pole separation give a finite compact-strip entire-multiplier bound of polynomial degree nine and exponential type F(t)-F(s), F=atan(u)/2+u/(2*(1+u^2)). A fresh global positive tempered Gaussian scalar state has the exact state-independent causal Kubo commutator. No global continuation of the off-clock nearby solution, Hadamard stress, interacting control or UV completion is asserted.",
        "not_established":["Global continuation of the separately certified off-clock nearby solution or nonlinear stability of either background",
            "Canonical interaction, omitted-operator, loop, threshold or quantitative interacting validity-band bounds",
            "Hadamard state, renormalized stress, semiclassical backreaction, interacting microcausality or global Fock unitary implementation",
            "Physical vacuum, amplitudes, finite-gravity dispersion and actual common-parent V/G/B matching",
            "Original realistic-field P8(a), original UV classification P8(b) or original P8 closure"],
        "verification_boundary":"Native literal constrained Hamiltonian and source, exact degree-two polynomial-q density generator, complete two-chart Euler and 4x4 phase intertwining identities, actual global gradient/kinetic bridge, rational large-complex-frequency degree and denominator-factor checks, explicit chart overlap and pole separation, exact real-frequency Hermitian cancellation, Gaussian Gram/CCR and Kubo identities. Written compact-strip logarithmic-norm composition, original entire continuation, distributional support and global Schwartz/Duhamel state arguments are pinned. Prefactors are proved finite on each compact time strip, not numerically bounded uniformly for all time. Every ancestor, source/proof/test byte and report field is independently replayed without scientific-library patches. Not proof-assistant formalization."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The global retuned scalar-support report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.98.GLOBAL_RETUNED_SCALAR_SUPPORT replay passed; interacting quantum V/G/B and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
