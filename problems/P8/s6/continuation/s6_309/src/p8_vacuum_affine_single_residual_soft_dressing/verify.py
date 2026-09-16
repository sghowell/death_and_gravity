"""Read-only same-state leading-soft dressing of one finite residual."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_leading_ladder_coulomb_phase import verify as previous
from p8_vacuum_affine_radiative_angular_finite import verify as radiative_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-single-residual-soft-dressing.json"
)
PARENT_SHA = "b9698ae4dbd703253e7d6fdd011d00370fed6c72d256bfb246bc6a0ee216f1cf"
RADIATIVE_SHA = "5e18135d1a49c2846047f2a71c3de00f7c82b84a131d8a6f924890308f4ccd99"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_single_residual_soft_dressing/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen full-D leading-ladder parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(radiative_input.REPORT) != RADIATIVE_SHA:
        raise ValueError("The frozen finite radiative angular conversion changed")
    radiative_input.validate_report(
        json.loads(radiative_input.REPORT.read_text()), radiative_input.build_report()
    )
    return {
        "S6_308_complete_full_D_leading_ladder_phase_rebuilt": PARENT_SHA,
        "S6_301_complete_finite_radiative_angular_conversion_rebuilt": RADIATIVE_SHA,
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
        raise ValueError("A single-residual soft-dressing proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.309.STATE_CORRECT_ALL_N_LEADING_SOFT_DRESSING_OF_ONE_FINITE_RESIDUAL",
        "date": "2026-09-16",
        "status": "SCOPED_EXACT_SINGLE_RESIDUAL_LEADING_SOFT_DRESSING; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/series.md",
            "notes/marked.md",
            "notes/physical.md",
            "notes/comparison.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_signed_soft_series": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_marked_operator_and_weighted_variation_bound": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged full47-tree signed one-graviton remainder has a complete same-state all-N additional-leading-soft calorimetric dressing. A uniform absolute bound on the signed dimensional series justifies its regulator limit and remaining-energy convolution. State-correct comparison with the unexpanded elastic reference bounds its relative magnitude by2e-768 and by a uniformly vanishing linear/quadratic threshold estimate. This named single-residual reference is not the full interacting multi-real rate; physical matching, further nonleading amplitudes, Regge and original V/G/B/P8 remain OPEN.",
        "not_established": [
            "The full all-N nonleading physical detector rate or an exact residual decomposition",
            "Two-soft contact terms, multiple finite residuals and correlated exact recoil",
            "Finite radiative hard loops, evanescent hard terms and physical matching",
            "A positive detector measure, unitarity or high-energy complex Regge control",
            "The original quantum state, common-parent bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact Bose and energy-simplex identities, signed-series uniform domination, separate positive-index and zero-index limit proofs, original47-tree signed-kernel calibrations, and weighted total-variation comparison. The detector power and remaining-energy suppression remain unexpanded. Pointwise positivity of a named reference is not probability-measure or full-QFT closure. Native/direct/ordinary/CLI use original SymPy; the adapter is FULL-only. Frozen ancestors and scoped P8(a) remain unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The single-residual soft-dressing report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.309 single-residual soft-dressing replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
