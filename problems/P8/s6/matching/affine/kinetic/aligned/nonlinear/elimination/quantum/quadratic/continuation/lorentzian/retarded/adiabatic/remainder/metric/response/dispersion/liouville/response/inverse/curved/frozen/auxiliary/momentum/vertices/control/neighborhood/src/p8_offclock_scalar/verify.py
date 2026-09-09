"""Read-only same-action off-clock coupled local-symbol replay."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_coupled_energy import verify as parent

from . import audit, datum, euler, principal, spatial

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"offclock-coupled-symbol.json"
PARENT_SHA="d00ab1e71e1a5a2e747d9efb00924c2b00fbbda9cb60b2566cdb1ba0fcdc51b6"
sha,serialize=affine.sha,affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_offclock_scalar/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen seven-mode hard-tree report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_78_fully_rebuilt":PARENT_SHA,
            "same_classical_action_actual_constraints_and_on_clock_Euler_control_replayed":True}


@cache
def build_report():
    prior=prior_checks()
    identities,gates=audit.residuals(),audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("An off-clock coupled-symbol gate failed")
    exact=affine.certify_residuals(identities)
    return {"schema":1,"claim":"P8-S6.79.OFFCLOCK_COUPLED_SYMBOL","date":"2026-09-08",
            "status":"CONSTRAINT_COMPATIBLE_NEARBY_DATA_AND_EULER_FIRST_HIGH_FREQUENCY_FROZEN_OBSTRUCTION; ACTUAL_INSTABILITY_CUTOFF_AND_ORIGINAL_P8_OPEN",
            "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
            "formulation":"FORMULATION.md",
            "written_proofs":["notes/datum.md","notes/spatial.md","notes/high-frequency.md","notes/scope.md"],
            "exact_residuals":exact,"named_exact_check_count":len(identities),
            "checked_scalar_entries":sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in identities.values()),
            "proof_checks":gates,
            "actual_constraint_compatible_homogeneous_family":serialize(datum.data()),
            "continuous_offclock_auxiliary_domain":serialize(datum.domain()),
            "complete_six_channel_phase_Hessian":serialize(spatial.phase_matrix()),
            "actual_canonical_gamma_blocks":serialize(principal.blocks()),
            "actual_Lagrangian_high_frequency_orders":serialize(euler.actual_orders()),
            "arbitrary_finite_time_jet_determinant_check":serialize(euler.generic_determinant()),
            "Euler_first_frozen_diagnostic":serialize(euler.result()),
            "old_clock_negative_controls":serialize({
                "incorrect_freeze_before_Euler":euler.result()["clock_negative_control_without_Euler_time_jets"],
                "actual_Euler_first_clock":euler.result()["actual_clock_Euler_leading_determinant"]}),
            "canonical_and_frame_conventions":"On the central flat hat slice, physical h=N*hat_h and physical squared momentum is q/N. W_i=partial_i sigma, P_sigma=-div(Pi), v=P_b/(2q), p=-2q*b. Keep the gamma time generator and actual finite homogeneous time jets before freezing. The physical matter coupling is not changed.",
            "quantifiers_and_boundary":"For every fixed real 0<|N-1|<=10^-28, there is a sufficiently large q beyond which the Euler-first frozen three-scalar constant determinant is negative and the kinetic matrix positive. Continuity in positive real growth frequency gives a frozen growing root. No quantitative q threshold, actual time-dependent instability, controlled interacting cutoff or UV exclusion is established.",
            "controls":audit.controls(),
            "verdict":"The original action admits an explicit constraint-compatible nearby homogeneous family inside its quantified auxiliary domain. Complete spatial reduction retains a scalar-momentum/longitudinal-vector coupling absent exactly on the clock. Its Legendre square supplies a negative q^4 Euler determinant coefficient that arbitrary finite background time jets cannot cancel. The original clock Euler sign remains positive and a false freeze-before-variation control is detected. This is an off-clock high-frequency frozen-neighborhood obstruction, not instability of the actual bounce or a UV verdict.",
            "not_established":["Numerical frequency threshold in a controlled EFT band or an interacting Wilsonian cutoff",
                "Actual time-dependent or nonlinear spatial instability from the frozen positive-real root",
                "Omitted higher-operator/loop/threshold matching errors or a quantum state theorem off the clock",
                "Transfer or contradiction of the exact on-clock finite-tree and separate quantum-response examples",
                "Finite Wilson matching, common-parent V/G/B, UV exclusion/completion or original P8 closure"],
            "verification_boundary":"Exact same-action lapse, spatial constraint and phase identities; continuous rational enclosures; written homogeneous smoothness, momentum-order and Euler-first determinant/continuity proofs. Not proof-assistant formalization."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The off-clock coupled-symbol report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.79.OFFCLOCK_COUPLED_SYMBOL replay passed; actual instability, cutoff and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
