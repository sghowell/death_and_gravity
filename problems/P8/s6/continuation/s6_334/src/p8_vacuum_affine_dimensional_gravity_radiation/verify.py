"""Read-only dimensionally complete selected radiation-tree report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_joint_soft_regulator import verify as input0
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import verify as input2
from p8_vacuum_affine_minimal_gravity_radiation import verify as input1
from p8_vacuum_affine_physical_virtual_soft_pairing import verify as input3

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-dimensional-gravity-radiation.json"
)
PARENT_SHA = "dfd88eebf2e1be8597025dc4f712e2b3fd8d86d5fac13de472c78271f6d82fc2"
GRAVITY_SHA = "a34d5160c55891304eda5f83bf90970854774e8aae327d9375b89e5fa797063e"
MATTER_SHA = "c5aeb643e173157f280201af2276924d5a211a5f1f97a271a0ed8dfc8c6d4d5f"
PHASE_SHA = "40c04cadfa1a1c542c16eed48f6df12bcfe0605ad69b39dbb70fcefaefe6b3ed"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_dimensional_gravity_radiation/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    for module, digest in (
        (input0, PARENT_SHA),
        (input1, GRAVITY_SHA),
        (input2, MATTER_SHA),
        (input3, PHASE_SHA),
    ):
        if sha(module.REPORT) != digest:
            raise ValueError("A frozen dimensional-tree input changed")
        module.validate_report(
            json.loads(module.REPORT.read_text()), module.build_report()
        )
    return {
        "S6_333_quantitative_joint_additional_soft_regulator_boundary_rebuilt": PARENT_SHA,
        "S6_304_canonical21_gravity_vertices_and_general_four_D_Ward_rebuilt": GRAVITY_SHA,
        "S6_295_existing_general_D_complete26_matter_tensor_rebuilt": MATTER_SHA,
        "S6_296_existing_rank_two_sew_and_dimensional_phase_boundary_rebuilt": PHASE_SHA,
        "S279_remains_rejected_and_archived": True,
        "S275_S276_S277_and_scoped_P8a_unchanged": True,
        "all6_historical_physical_qualifications_unchanged": True,
        "all_original_primitive_and_matching_frontiers_retained": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A dimensional radiation proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.334.DIMENSIONALLY_COMPLETE_GRAVITY_RADIATION_TREE_AND_POLARIZATION_SEW",
        "date": "2026-09-18",
        "status": "SCOPED_EXACT_DIMENSIONAL_TREE_AND_SEW; NOT_FINITE_HARD_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/vertices.md",
            "notes/decomposition.md",
            "notes/polarizations.md",
            "notes/bounds.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_dimensionally_complete_tree": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_canonical_47_graph_tensor_and_polarization_sew": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The original selected47 tree has a general-D Ward identity and exact gravity correction 2mu^2*(D-4)/(D-2)*S_aux. Its complete polarization sew is cubic in epsilon/(1+epsilon), retaining the trace contribution invisible to D4 helicities and all matter-gravity interference. A pointwise explicit regulator-rate majorant is established. Auxiliary exchange is an algebraic identity only. Integrated gravity real-minus-soft pairing, finite hard-loop matching and original P8 closure are not established.",
        "not_established": [
            "Integrated full-gravity real-minus-soft dimensional limit and virtual hard pairing",
            "Finite hard-loop evanescent matching or loop-corrected radiation",
            "Full interacting detector probability or all-N hard summability",
            "Quantum unitarity and absolute complex Regge control",
            "Common-parent bounce, UV or original V/G/B/P8 closure",
        ],
        "verification_boundary": "General invariant-index rational Ward/decomposition identities and written null-rank-one span and polarization proofs. Independent integer-D canonical component expansions and literal matter graphs calibrate every polarization, with rational fixtures and a separate unchanged-original-point check. Finite checks do not replace the symbolic proof. Not kernel-formalized; original SymPy outside FULL.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The dimensional radiation report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.334 complete dimensional radiation-tree replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
