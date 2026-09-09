"""Read-only actual seven-mode column and finite hard-tree replay."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_coupled_vertices import verify as parent

from . import audit, bounds, energy, other_modes, scalars, tree, window

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"seven-mode-hard-tree.json"
PARENT_SHA="37ebcea32c48637d390faeacc51c5e2074ab4c380e6fbb504f25f2674a446f8c"
sha,serialize=affine.sha,affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_coupled_energy/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen coupled higher-vertex report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_77_fully_rebuilt":PARENT_SHA,
            "full_classical_seven_mode_H3_H4_and_original_selected_vector_state_replayed":True}


@cache
def build_report():
    prior=prior_checks()
    identities,gates=audit.residuals(),audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A seven-mode finite-tree gate failed")
    exact=affine.certify_residuals(identities)
    return {"schema":1,"claim":"P8-S6.78.SEVEN_MODE_HARD_TREE","date":"2026-09-08",
            "status":"ACTUAL_SEVEN_MODE_FREE_COLUMN_AND_FINITE_HARD_TREE_BOUNDS; INTERACTING_WILSONIAN_CUTOFF_AND_ORIGINAL_P8_OPEN",
            "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
            "formulation":"FORMULATION.md",
            "written_proofs":["notes/scalar-energy.md","notes/columns.md","notes/other-modes.md","notes/tree.md"],
            "exact_residuals":exact,"named_exact_check_count":len(identities),
            "checked_scalar_entries":sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in identities.values()),
            "proof_checks":gates,
            "actual_two_chart_scalar_matrices":serialize({name:scalars.data(name) for name in ("unitary","gamma")}),
            "actual_background_jet_enclosures":serialize(bounds.backgrounds()),
            "continuous_chart_margins":serialize(bounds.elementary_margins()),
            "canonical_initial_shift_negative_control":serialize(energy.initial_mixing_negative_control()),
            "normalized_scalar_column_domain":serialize(window.data()),
            "unchanged_selected_vector_columns":serialize(other_modes.vectors()),
            "both_tensor_columns":serialize(other_modes.tensors()),
            "finite_hard_tree_estimate":serialize(tree.data()),
            "exact_phase_homothety":"With L=M*tau, every raw canonical coordinate AND momentum scales by 1/L; L^2 H_n(Q/L,p/L)=L^(2-n) H_n(Q,p). A linear time-dependent canonical map adds only a quadratic generator. The Hamiltonian's 1/tau cancels dt=tau*du.",
            "fixed_fiber_observable":"Finite-time 1-to-2 and 2-to-1 cubic blocks, and the hard-masked connected 2-to-2 tree through H4 and H3 squared, with d^3k/(2*pi)^3 on fixed-total-momentum fibers of R^3. All seven modes, both phase components, the fixed-center volume, both time orders and internal orientations are retained. Finite-band mode-kernel construction, not globally implementable homogeneous Fock evolution.",
            "controls":audit.controls(),
            "verdict":"The actual coupled scalar energy, tensor energy and unchanged selected-state vector comparison bound all seven canonical phase columns on the explicit short-window high-frequency domain. Full reduced H3/H4, hard transfers, Wick/species counts and a fixed-fiber Schur estimate yield exact C3/L and C4/L^2 bounds. The separately named L=10^353 example makes both <=1/1000. This is a finite-order classical hard-tree statement, not a Wilsonian cutoff, full Fock norm, UV completion or original P8 closure.",
            "not_established":["Uniform interacting Wilsonian/all-orders cutoff, higher loops/operators/threshold errors or nonlinear spatial existence",
                "Zero/forward channels, infinite-time scattering, a full Fock-space operator norm or globally implementable homogeneous evolution",
                "All-momentum scalar/tensor Hadamard states or new quantum stability/cones",
                "Transfer of the separate L=10^24 response or frozen-pole verdicts to this named hierarchy",
                "Finite Wilson matching, common-parent V/G/B or original P8 closure"],
            "verification_boundary":"Exact symbolic two-component identities and continuous rational chart bounds, with written canonical normalization, unchanged-state mode comparison and finite-tree kernel/Schur proofs. Not proof-assistant formalization."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The seven-mode finite hard-tree report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.78.SEVEN_MODE_HARD_TREE replay passed; Wilsonian cutoff and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
