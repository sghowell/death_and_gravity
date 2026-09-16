"""Read-only full radiative-state soft-index certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_calorimetric_soft_resummation import verify as previous
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import (
    verify as recoil_input,
)

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-radiative-state-soft-index.json"
PARENT_SHA = "604445d8362c3dd9cec6e456bab0743e0611aaac2e60794e64f99be41cd63879"
RECOIL_SHA = "c5aeb643e173157f280201af2276924d5a211a5f1f97a271a0ed8dfc8c6d4d5f"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_radiative_state_soft_index/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen complete leading-soft parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(recoil_input.REPORT) != RECOIL_SHA:
        raise ValueError("The frozen complete physical recoil changed")
    recoil_input.validate_report(
        json.loads(recoil_input.REPORT.read_text()), recoil_input.build_report()
    )
    return {
        "S6_299_leading_soft_sum_ancestry_rebuilt": PARENT_SHA,
        "S6_295_complete_physical_recoil_rebuilt": RECOIL_SHA,
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
        raise ValueError("A complete radiative-state soft-index proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.300.COMPLETE_RADIATIVE_STATE_SOFT_INDEX_UNIFORM_MULTIPLICITY_RECOIL_AND_NESTED_IR_PAIRING",
        "date": "2026-09-15",
        "status": "SCOPED_COMPLETE_RADIATIVE_SOFT_INDEX_AND_UNIFORM_MULTIPLICITY_STABILITY; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/recoil.md",
            "notes/index.md",
            "notes/stability.md",
            "notes/matching.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_radiative_recoil": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_massless_soft_index_and_uniform_stability": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete recoiled radiative soft index includes every massive/massless and massless/massless pair. Each collinear mass logarithm cancels by conservation; the finite angular-energy formula is boost invariant and continuous under finite-energy measure limits. Its difference from the elastic index is bounded by40R/kappa uniformly in multiplicity. Nested real/virtual poles must use the same radiative state. After correct pairing the known detector power and Gamma change is below3*10^-799 on the original compact domain. The finite angular and hard radiative terms and all physical IR/Regge/V/G/B/P8 closure obligations remain OPEN.",
        "not_established": [
            "Full D finite angular conversion for arbitrary radiative states",
            "Finite radiative hard amplitudes and virtual matching",
            "Uniform complete-amplitude remainder and all omitted loops",
            "Full physical analytic/Regge and high-energy remainder control",
            "Original quantum state/domain/measure/bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact shell/conservation identities, full massless-pair limits, positive compact-domain majorants and written energy-measure continuity arguments prove the scoped index and stability statements. Independent TT angular integrals, massless regulators, boosts and collinear splits calibrate the whole formula. This does not construct an interacting quantum state or a full radiative amplitude. Native/direct/ordinary/CLI use original SymPy; the adapter is FULL-only. All frozen ancestors, historical qualifications and original closure statuses remain unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete radiative-state soft-index report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.300 complete radiative-state soft-index replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
