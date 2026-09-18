"""Read-only uniform radiative logarithmic-coefficient report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_minimal_gravity_radiation import verify as action_input
from p8_vacuum_affine_radiative_soft_state_transfer import verify as previous
from p8_vacuum_affine_radiative_state_soft_index import verify as state_input
from p8_vacuum_affine_uniform_all_tree_bound import verify as tree_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT
    / "certificates/polynomial-vacuum-affine-radiative-logarithmic-soft-coefficient.json"
)
PARENT_SHA = "32396aad07a5c1521a282fe6161c0a8e99ff6e312f55df4937de8f5501e30655"
STATE_SHA = "83fc98b6541ae6fb75c18c342d20c3e4bc03842fe69a3191b3dc8d1a09509df7"
ACTION_SHA = "a34d5160c55891304eda5f83bf90970854774e8aae327d9375b89e5fa797063e"
TREE_SHA = "6239bb0b855c48275f8598cdfdec3f8dfa2590545d5b4abb5b6ddfe765017b0b"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(
            ROOT.glob(
                "src/p8_vacuum_affine_radiative_logarithmic_soft_coefficient/*.py"
            )
        )
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    for module, digest in (
        (previous, PARENT_SHA),
        (state_input, STATE_SHA),
        (action_input, ACTION_SHA),
        (tree_input, TREE_SHA),
    ):
        if sha(module.REPORT) != digest:
            raise ValueError("A frozen radiative logarithmic input changed")
        module.validate_report(
            json.loads(module.REPORT.read_text()), module.build_report()
        )
    return {
        "S6_325_uniform_radiative_soft_state_continuity_rebuilt": PARENT_SHA,
        "S6_300_exact_radiative_recoil_and_finite_index_rebuilt": STATE_SHA,
        "S6_304_original_Einstein_scalar_action_and_full_N1_rebuilt": ACTION_SHA,
        "S6_319_complete_all_finite_tree_majorant_rebuilt": TREE_SHA,
        "S279_remains_rejected_and_archived": True,
        "S275_S276_S277_and_scoped_P8a_unchanged": True,
        "all6_historical_physical_qualifications_unchanged": True,
        "all_original_primitive_and_matching_frontiers_retained": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A radiative logarithmic-coefficient proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.326.UNIFORM_RADIATIVE_LOGARITHMIC_SOFT_COEFFICIENT",
        "date": "2026-09-17",
        "status": "SCOPED_UNIFORM_RADIATIVE_LOG_COEFFICIENT; NOT_FULL_LOOP_ALL_N_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/massive.md",
            "notes/conservation.md",
            "notes/null.md",
            "notes/bounds.md",
            "notes/measure.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_logarithmic_conventions": serialize(
            {name: payload(packets[name]) for name in names[:3]}
        ),
        "whole_radiative_coefficient_bounds_and_named_measure": serialize(
            {name: payload(packets[name]) for name in names[3:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The attributed proved one-loop logarithmic soft coefficient has a collinear-finite expression on every finite original radiative state, with |C|<320 and |C-C_Born|<25000R. The physical coefficient scales as kappa^(-3/2). Differentiation of conservation-zero terms is retained until cancellation of auxiliary masses and individual energy logarithms with the phase. Exact collinear energy splitting leaves the result unchanged; exactly soft-collinear atoms receive no unique point value. The named logarithmic addition has integrable signed interference with the complete original47-tree amplitude. No full hard loop, uniform finite remainder, common regulator, all-N hard measure or original P8 closure is claimed.",
        "not_established": [
            "Finite hard loop and uniform O(1) one-loop soft remainder",
            "Detector-scale term and complete common-regulator real-virtual matching",
            "Simultaneous multi-soft remainder, all-N hard measure or positive normalized rate",
            "Interacting quantum state, unitarity or absolute complex Regge",
            "Original common-parent bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Original generic scalar vertex/prescription and Born-current convention audit; removable real kernel and retained imaginary phase; exact differentiation-before-conservation identities; original rational recoil and complex TT calibrations; exact finite collinear splitting and nonunique soft-collinear counterexample; analytic all-finite-multiplicity coefficient/difference bounds and selected signed logarithmic integrals. The written analytic argument supplies uniformity, not sample tests or arithmetic gates alone. The external theorem is attributed, not rederived from every loop diagram. Not kernel-formalized; original SymPy outside FULL.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The radiative logarithmic report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.326 uniform radiative logarithmic coefficient replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
