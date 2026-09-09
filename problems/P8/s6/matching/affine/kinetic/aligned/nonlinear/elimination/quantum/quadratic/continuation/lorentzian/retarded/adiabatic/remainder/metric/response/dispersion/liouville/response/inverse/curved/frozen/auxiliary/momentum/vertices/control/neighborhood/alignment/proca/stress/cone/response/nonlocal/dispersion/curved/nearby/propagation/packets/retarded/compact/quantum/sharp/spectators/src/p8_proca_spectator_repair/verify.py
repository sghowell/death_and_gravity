"""Read-only finite principal-extension and actual Proca spectator screening certificate."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as original
from p8_proca_sharp_response import verify as parent

from . import audit, pencil, proca, repair

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"finite-principal-spectator-repair.json"
PARENT_SHA="95f706f2efd91d4bdb418d32ab29451c63f73b1cbcf4969b0b339ff7e7892380"
sha,serialize=original.sha,original.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_proca_spectator_repair/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen sharper actual compact-response report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_94_fully_rebuilt":PARENT_SHA,
        "same_actual_nearby_light_action_background_and_physical_matter_frame":True,
        "ordinary_Proca_classical_application_retains_its_own_constraints":True,
        "no_effective_light_block_assumed_to_be_a_bare_parent_restriction":True}


@cache
def build_report():
    prior=prior_checks()
    identities,gates=audit.residuals(),audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A finite principal-extension or Proca screening gate failed")
    exact=original.certify_residuals(identities)
    return {"schema":1,"claim":"P8-S6.95.FINITE_PRINCIPAL_SPECTATOR_REPAIR","date":"2026-09-09",
        "status":"FINITE_POSITIVE_KINETIC_HERMITIAN_PRINCIPAL_EXTENSIONS_NEED_ACTUAL_LIGHT_K_H_REPAIR_IN_BOTH_ORIENTATIONS; ORDINARY_PROCA_SPECTATORS_RETAIN_FAST_SCALAR_FACTOR; GENERAL_UV_AND_ORIGINAL_P8_OPEN",
        "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
        "formulation":"FORMULATION.md",
        "written_proofs":["notes/principal.md","notes/repair.md","notes/proca.md","notes/scope.md"],
        "exact_residuals":exact,"named_exact_check_count":len(identities),
        "checked_scalar_entries":sum(value.rows*value.cols if isinstance(value,sp.MatrixBase) else 1
                                     for value in identities.values()),
        "proof_checks":gates,
        "actual_fixed_physical_matter_frame_light_canonical_forms":serialize(pencil.light()),
        "finite_Hermitian_principal_extension_and_explicit_embedding":serialize(pencil.extension()),
        "necessary_two_orientation_light_repair_and_negative_controls":serialize(repair.data()),
        "actual_ordinary_Proca_constraint_decoupling_and_finite_species":serialize(proca.data()),
        "controls":audit.controls(),
        "verdict":"The actual nearby light principal pair is I and diag(c_clock^2,1) in its old-kinetic canonical comparison, with c_clock^2-1>5*10^-6 throughout the interval. For any finite Hermitian second-order physical principal pencil with positive full kinetic form and this explicit light embedding, unchanged light kinetic and gradient blocks leave a real characteristic of absolute speed at least c_clock, including arbitrary finite mixed time-space and heavy/cross blocks. Matter-causal propagation in both orientations requires the sum of relative light kinetic and gradient change norms to exceed 5*10^-6. This is only necessary. The literal ordinary-Proca addition at zero vector background decouples quadratically, and its temporal constraint gives three healthy physical modes with luminal principal cone. Any finite number of such classical spectators leaves the faster scalar factor. Genuine light-action changes and frequency-dependent parent response remain separate research routes; no general UV no-go is proved.",
        "not_established":["A covariant light-action repair, new on-shell solution or healthy interacting completion",
            "An explicit bare-parent light embedding or positive physical kinetic form for an arbitrary UV theory",
            "Higher-derivative, nonlocal, loop-induced, infinite-tower or unreduced constrained extension classification",
            "Joint temporal/spatial interacting EFT validity, actual parent matching errors or that parent's own response tail",
            "Hadamard stress, semiclassical backreaction, V/G/B, realistic-field P8(a) or original P8 closure"],
        "verification_boundary":"Native exact fixed-frame canonical congruence, complete mixed principal pencil restriction, both orientations, literal vector/metric quadratic variation, constant canonical mass normalization, temporal Proca constraint and finite lower-order mass anchor. Written Hermitian continuity and operator-norm proofs cover arbitrary finite dimension; written decoupling proves the all-species product. Exact domain controls reject invalid budgets and species. Every parent is natively rebuilt and every source/proof/test byte pinned. No scientific-library patch, inferred bare-parent embedding, assigned UV error or frozen action mutation. Not proof-assistant formalization."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The finite principal spectator-repair report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.95.FINITE_PRINCIPAL_SPECTATOR_REPAIR replay passed; changed light action, general UV and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
