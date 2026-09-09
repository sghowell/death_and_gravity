"""Read-only isolated ordinary-Proca rank-one massive scalar inverse."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as original
from p8_proca_nonlocal_response import verify as parent

from . import audit, cut, kernel, spectral

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"rank-one-massive-inverse.json"
PARENT_SHA="fea66d1775ce50b7ef941af5b0987269f86ab6f1e6a5a3a23fc53ca63ba3bb01"
sha,serialize=original.sha,original.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_proca_rank_one_inverse/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen prepared ordinary Proca response report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_85_fully_rebuilt":PARENT_SHA,
            "new_rank_one_reference_not_old_full_matrix_inverse_transferred":True}


@cache
def build_report():
    prior=prior_checks()
    identities,gates=audit.residuals(),audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A new rank-one massive scalar inverse gate failed")
    exact=original.certify_residuals(identities)
    return {"schema":1,"claim":"P8-S6.86.RANK_ONE_MASSIVE_INVERSE","date":"2026-09-08",
            "status":"NEW_ISOLATED_RANK_ONE_MASSIVE_RANGE_FACTOR_ZERO_FREE_WITH_POSITIVE_CUT_AND_L1_CAUSAL_SCALAR_INVERSE; ACTUAL_COUPLED_INVERSE_MATCHING_AND_ORIGINAL_P8_OPEN",
            "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
            "formulation":"FORMULATION.md",
            "written_proofs":["notes/block.md","notes/first-sheet.md","notes/cut.md","notes/kernel.md","notes/scope.md"],
            "exact_residuals":exact,"named_exact_check_count":len(identities),
            "checked_scalar_entries":sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in identities.values()),
            "proof_checks":gates,
            "actual_new_physical_flat_pair_and_fixed_canonical_contacts":serialize({
                sector:{"pair":spectral.pair(sector),"contact":spectral.contact(sector)} for sector in ("T","L")}),
            "new_exact_rank_one_massive_scalar_block_and_radial_representation":serialize(spectral.data()),
            "full_massive_positive_axis_and_removable_zero_controls":serialize({
                value:spectral.positive_axis(value) for value in (0,1,100)}),
            "new_positive_scalar_inverse_cut_density_and_exact_moments":serialize(cut.data()),
            "exact_open_cut_controls":serialize({str(value):cut.fraction(value) for value in
                (sp.Rational(1,10),sp.Rational(1,2),sp.Rational(9,10))}),
            "new_causal_scalar_L1_kernel_and_primitive":serialize(kernel.data()),
            "positive_fixed_mass_moment_and_primitive_bounds":serialize({value:kernel.moment_bound(value) for value in (1,2,1000)}),
            "controls":audit.controls(),
            "verdict":"The new physical reference is rank one, not invertible as a two-source matrix. Its exact scalar factor F=-H has Re(H)>=16/15 on the entire first sheet, with no threshold, cut-bank or finite first-sheet zero. Its inverse has a positive cut density with exact moments 1/4 and 9/(560 m^2), no instantaneous or pole remainder, and an ordinary L1 retarded time kernel. The complete massive formula excludes the spurious pole suggested by a high-frequency logarithmic truncation. Only the normalized scalar range factor is inverted.",
            "not_established":["A full two-source inverse, actual curved selected-state or coupled tree-plus-loop inverse",
                "A stationary healthy Minkowski vacuum of the full candidate or physical propagator spectral positivity",
                "Quantum stability/cones, a nonlinear neighborhood or independently varied initial-state response",
                "Canonical profile interactions, mixed loops, omitted operators, heavy thresholds or Wilsonian cutoff",
                "Finite Wilson matching, common-parent V/G/B, UV completion/exclusion or original P8 closure"],
            "verification_boundary":"Exact new canonical contact, adiabatic Taylor, spectral moment, threshold, imaginary-sign, cut and normalization identities with written harmonic minimum-principle, Cauchy-contour and oscillatory-kernel proofs. The L1 existence proof does not invent a numerical L1 norm. Every parent is natively rebuilt. Not proof-assistant formalization."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The rank-one massive inverse report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.86.RANK_ONE_MASSIVE_INVERSE replay passed; actual coupled inverse, matching and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
