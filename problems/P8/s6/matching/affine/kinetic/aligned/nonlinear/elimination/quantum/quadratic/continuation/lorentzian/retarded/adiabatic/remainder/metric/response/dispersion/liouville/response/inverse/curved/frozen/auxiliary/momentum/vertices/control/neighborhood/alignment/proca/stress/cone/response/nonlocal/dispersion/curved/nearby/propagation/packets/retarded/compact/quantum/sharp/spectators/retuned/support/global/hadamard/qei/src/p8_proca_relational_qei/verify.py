"""Read-only actual relational quadratic bounds and physical-stress boundary."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as certificate
from p8_proca_global_hadamard import verify as parent

from . import audit, modes, observers, qei, stress

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"relational-quadratic-inequalities.json"
PARENT_SHA="e194ade04ac78efa48f1f73165c07a2bf2ce1f9eef234a66aa72d89ebb55dd41"
sha,serialize=certificate.sha,certificate.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_proca_relational_qei/*.py"))+sorted(ROOT.glob("tests/*.py"))
        +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen global two-cone Hadamard-state report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_99_fully_rebuilt":PARENT_SHA,
        "same_actual_global_action_background_and_generalized_state_class":True,
        "no_free_matter_or_electromagnetic_state_substituted":True}


@cache
def build_report():
    prior=prior_checks()
    identities,gates=audit.residuals(),audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A relational mode, observer or quadratic-bound gate failed")
    exact=certificate.certify_residuals(identities)
    return {"schema":1,"claim":"P8-S6.100.RELATIONAL_QUADRATIC_INEQUALITIES","date":"2026-09-09",
        "status":"COMPLETE_ACTUAL_GLOBAL_RELATIONAL_SCALAR_SUBCLOCK_REFERENCE_DIFFERENCE_QUADRATIC_BOUNDS_AND_LEADING_SHORT_SAMPLING_COEFFICIENTS_WITH_NONZERO_CLOCK_RESIDUE_AND_EXPLICIT_TEST_STRESS_CONSERVATION_OBSTRUCTION; COVARIANT_PHYSICAL_STRESS_COSMOLOGICAL_REMAINDERS_V_G_B_AND_ORIGINAL_P8_OPEN",
        "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
        "formulation":"FORMULATION.md",
        "written_proofs":["notes/modes.md","notes/observers.md","notes/inequalities.md","notes/short-sampling.md","notes/stress.md","notes/scope.md"],
        "exact_residuals":exact,"named_exact_check_count":len(identities),
        "checked_scalar_entries":sum(value.rows*value.cols if isinstance(value,sp.MatrixBase) else 1
                                     for value in identities.values()),
        "proof_checks":gates,
        "actual_normalized_modes_and_nonzero_relational_clock_residue":serialize(modes.data()),
        "actual_worldline_domain_and_interluminal_control":serialize(observers.data()),
        "finite_reference_difference_QEI_and_short_sampling_coefficients":serialize(qei.data()),
        "fixed_background_test_stress_and_actual_coupled_source":serialize(stress.data()),
        "controls":audit.controls(),
        "verdict":"The actual reduced relational scalar admits finite reference-state difference inequalities for real differential squares on every compactly sampled worldline slower than the clock mode; the all-time physical-speed class |v|<=.99 is explicit. The exact normalized clock spectral weight is nonzero, so a free single-matter-metric Hadamard replacement cannot reproduce the actual commutator. At the comoving bounce the leading short-sampling test-energy coefficient is hbar/(16*pi^2*kappa) times 1+10863*sqrt(17985)/1723683599, between 1.00084 and 1.00085 times the unmixed reference coefficient. Full finite-width constants remain reference-state functionals. The fixed-background scalar test tensor has an actual nonzero coupled clock divergence source and is not declared the complete conserved physical stress. No optimal QEI, interluminal no-go, cosmological remainder or original P8 closure follows.",
        "not_established":["Optimal or saturated QEI, numerical finite-width reference functional, or controlled cosmological-scale short-sampling remainder",
            "Nonexistence of every interluminal worldline restriction or a state-independent interluminal QEI no-go",
            "Conserved full covariant renormalized physical stress, interacting Ward identities, semiclassical backreaction or realistic-field singularity theorem",
            "Full quantum tensor/vector system, interacting/omitted-operator/loop/threshold validity or actual vacuum/dispersion/common-parent V/G/B matching",
            "Original null/nonminimal/realistic-field P8(a), original UV classification P8(b) or original P8 closure"],
        "verification_boundary":"Native normalized two-mode symplectic/kinetic bridge, complete order-zero transport and all basis time jets, nonzero actual matter clock residue, whole-interval physical observer contraction and interluminal controls, exact radial frequency moments and algebraic short-sampling coefficients, arbitrary-gradient/Hessian covariant test-stress divergence and actual coupled scalar principal source. Pinned written elliptic clock-wavefront, positive-type pullback, reference difference inequality, asymptotic sampling-limit and covariant-source scope proofs. All ancestors, source/proof/test bytes and report fields are natively replayed. Scientific SymPy is unchanged. Written distribution and limiting arguments are not proof-assistant formalization."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The relational quadratic-inequality report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.100.RELATIONAL_QUADRATIC_INEQUALITIES replay passed; physical stress, cosmological remainders, V/G/B and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
