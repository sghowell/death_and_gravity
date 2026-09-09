"""Read-only sharper actual nearby transfer and rational compact-response certificate."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as original
from p8_proca_scalar_ccr import verify as parent

from . import audit, basis, jets, probes, transfer

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"sharper-rational-compact-response.json"
PARENT_SHA="405fadc4468a18123330a9607d7c5117aea24db0df1f7a5818e7dfd65e5ba0b3"
sha,serialize=original.sha,original.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_proca_sharp_response/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen actual nearby quadratic scalar CCR report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_93_fully_rebuilt":PARENT_SHA,
        "same_actual_background_full_quadratic_action_retarded_kernel_and_CCR_algebra":True,
        "new_named_probe_pair_not_an_overwrite_of_old_widths":True,
        "no_interacting_EFT_cutoff_or_parent_matching_transferred":True}


@cache
def build_report():
    prior=prior_checks()
    identities,gates=audit.residuals(),audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A sharper actual complex-jet or response gate failed")
    exact=original.certify_residuals(identities)
    return {"schema":1,"claim":"P8-S6.94.SHARPER_RATIONAL_COMPACT_RESPONSE","date":"2026-09-09",
        "status":"ACTUAL_COMPLEX_JETS_GIVE_SHARPER_COMPLETE_TRANSFER_FRONT_REMAINDER_AND_NEW_RATIONAL_COMPACT_PROBES_WITH_SPATIAL_CUTOFF_BELOW_1E94_AND_NONZERO_QUADRATIC_COMMUTATOR; TEMPORAL_EFT_INTERACTING_PARENT_MATCHING_AND_ORIGINAL_P8_OPEN",
        "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
        "formulation":"FORMULATION.md",
        "written_proofs":["notes/jets.md","notes/basis.md","notes/transfer.md","notes/probes.md","notes/scope.md"],
        "exact_residuals":exact,"named_exact_check_count":len(identities),
        "checked_scalar_entries":sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in identities.values()),
        "proof_checks":gates,
        "actual_complex_whole_box_logarithmic_jets":serialize(jets.data()),
        "small_exact_action_basis_and_complete_source_output_factor":serialize(basis.data()),
        "complete_literal_finite_frequency_Laurent_majorants":serialize(transfer.coefficient_majorants()),
        "sharper_all_frequency_transfer_and_two_endpoint_normal_form":serialize(transfer.data()),
        "new_rational_compact_probes_complete_remainder_and_spatial_band":serialize(probes.data()),
        "controls":audit.controls(),
        "verdict":"On the unchanged actual nearby local bounce, direct complex-box logarithmic derivatives bound the complete action-normalized basis by 6, its inverse by 4 and its derivative by 1/50. Exact real-time packet transfer is below 48 for all k>=1; the unit ball uses the regular original polynomial generator. The same exact high-frequency normal form now gives a two-endpoint scalar front error below 10^16/k^2 for k>=10^14. Three complete frequency regions give spatial remainder L2 norm below 10^9. A newly named smooth compact probe pair with exact rational widths has full classical pairing at least 3D/4 and a spatial-band contribution at least D/2 below an explicit rational cutoff smaller than 10^94. Its compact matter-spacelike central commutator in the unchanged quadratic CCR algebra has magnitude at least 3*hbar*D/(4*kappa). No old probe, background or state is reassigned. This is not temporal/interacting EFT validity or actual UV-parent matching.",
        "not_established":["A joint temporal/spatial interacting EFT validity band, source apparatus or detector noise bound",
            "Actual common-parent low-band matching and that parent's own high-band response",
            "Canonical interactions, omitted operators, loops, heavy thresholds, finite Wilson matching or V/G/B",
            "Hadamard regularity, renormalized stress or nearby self-consistent semiclassical backreaction",
            "Nonperturbative local quantum-gravity observables, global nearby completeness, universal UV exclusion/completion or original P8 closure"],
        "verification_boundary":"Native exact whole-complex-box rational jet inequalities, literal compressed basis and full time derivative, complete Laurent matrices, Cauchy and gap-derivative constants, independent three-region Plancherel integrals, new rational compact widths, every support and plateau gate, and exact source/action-normalized quantum lower. Written holomorphic, evolution and compact-probe proofs are pinned. Every parent is natively rebuilt. No scientific-library patch, complex ordering, uncovered frequency interval, floating bound, old-probe mutation or cutoff/UV assumption. Not proof-assistant formalization."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The sharper rational compact-response report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.94.SHARPER_RATIONAL_COMPACT_RESPONSE replay passed; temporal EFT, interacting parent matching and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
