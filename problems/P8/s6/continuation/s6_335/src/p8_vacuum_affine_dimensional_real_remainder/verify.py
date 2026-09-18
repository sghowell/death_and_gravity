"""Read-only complete-tree dimensional real-minus-soft report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_dimensional_gravity_radiation import verify as input0
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import verify as input4
from p8_vacuum_affine_minimal_gravity_radiation import verify as input1
from p8_vacuum_affine_one_newton_inclusive_assembly import verify as input3
from p8_vacuum_affine_physical_virtual_soft_pairing import verify as input2

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-dimensional-real-remainder.json"
PARENT_SHA = "ce0ed27669fe9fefeb0064deae3a8565c2790ce517dae44f0737077f2dfc810b"
GRAVITY_SHA = "a34d5160c55891304eda5f83bf90970854774e8aae327d9375b89e5fa797063e"
PHASE_SHA = "40c04cadfa1a1c542c16eed48f6df12bcfe0605ad69b39dbb70fcefaefe6b3ed"
BORN_SHA = "6206f5d56959f01a8467e2c2b74ff43a2fad3cc75df0cee67e9bb8ce82f39564"
MATTER_SHA = "c5aeb643e173157f280201af2276924d5a211a5f1f97a271a0ed8dfc8c6d4d5f"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_dimensional_real_remainder/*.py"))
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
        (input2, PHASE_SHA),
        (input3, BORN_SHA),
        (input4, MATTER_SHA),
    ):
        if sha(module.REPORT) != digest:
            raise ValueError("A frozen dimensional real-remainder input changed")
        module.validate_report(
            json.loads(module.REPORT.read_text()), module.build_report()
        )
    return {
        "S6_334_complete_D_tree_and_canonical_trace_response_rebuilt": PARENT_SHA,
        "S6_304_complete_D4_tree_and_uniform_physical_recoil_gap_rebuilt": GRAVITY_SHA,
        "S6_296_actual_positive_fixed_domain_measure_and_matter_limit_rebuilt": PHASE_SHA,
        "S6_297_original_positive_full_Born_and_forward_boundary_rebuilt": BORN_SHA,
        "S6_295_original_recoil_and_canonical_matter_soft_remainder_rebuilt": MATTER_SHA,
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
        raise ValueError("A dimensional real-remainder proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.335.COMPLETE_TREE_DIMENSIONAL_REAL_MINUS_SOFT_LIMIT",
        "date": "2026-09-18",
        "status": "SCOPED_FULL_SELECTED_TREE_REAL_MINUS_SOFT_LIMIT; NOT_VIRTUAL_HARD_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/soft.md",
            "notes/embedding.md",
            "notes/bounds.md",
            "notes/continuity.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_same_D_Born_and_soft_remainder": serialize(
            {name: payload(packets[name]) for name in names[:3]}
        ),
        "whole_integrated_bounds_and_original_recoil_calibrations": serialize(
            {name: payload(packets[name]) for name in names[3:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete selected47 real-minus-soft integral tends to the S304 physical D4 remainder at fixed nonforward hard data and every0<x<=1/8, also jointly with removal of a lower energy cutoff. A six-dimensional trace completion and separate low/high-energy bounds retain all polarizations, original couplings and SAME-D Born subtraction. The fixed-resolution uniform bound is<[1e13+1e35*x^2/delta^2]/kappa, with a sharper vanishing low-cutoff bound. This is not a limit of separate divergent rates, a uniform forward theorem, finite virtual hard matching or original P8 closure.",
        "not_established": [
            "Uniform forward limit or IR finiteness of either unsubtracted real rate",
            "Finite virtual hard amplitude or hard-loop evanescent matching",
            "Full interacting detector probability or all-N hard summability",
            "Quantum unitarity and absolute complex Regge control",
            "Common-parent bounce, UV or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Written low/high-window estimates and fixed-domain dominated convergence use exact general soft/Born identities, D6 canonical trace completion and original recoil gaps. Ten exact original-source recoil configurations calibrate the formulas; finite sampling is not the uniform proof. All prior qualifications retained. Not kernel-formalized; original SymPy outside FULL.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The dimensional real-remainder report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.335 complete-tree dimensional real-remainder replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
