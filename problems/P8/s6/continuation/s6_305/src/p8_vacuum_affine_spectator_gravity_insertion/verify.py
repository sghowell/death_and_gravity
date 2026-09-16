"""Read-only complete fixed spectator insertion and known-rate certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_minimal_gravity_radiation import verify as previous
from p8_vacuum_affine_one_newton_inclusive_assembly import verify as born_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-spectator-gravity-insertion.json"
PARENT_SHA = "a34d5160c55891304eda5f83bf90970854774e8aae327d9375b89e5fa797063e"
BORN_SHA = "6206f5d56959f01a8467e2c2b74ff43a2fad3cc75df0cee67e9bb8ce82f39564"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_spectator_gravity_insertion/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen full Einstein radiation parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(born_input.REPORT) != BORN_SHA:
        raise ValueError("The frozen full positive Born normalization changed")
    born_input.validate_report(
        json.loads(born_input.REPORT.read_text()), born_input.build_report()
    )
    return {
        "S6_304_whole_Einstein_radiation_and_uniform_real_rate_rebuilt": PARENT_SHA,
        "S6_297_full_positive_Born_normalization_rebuilt": BORN_SHA,
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
        raise ValueError("A spectator-insertion proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.305.COMPLETE_FIXED_SPECTATOR_METRIC_INSERTIONS_AND_UNIFORM_SIGNED_KNOWN_REFERENCE_RATE",
        "date": "2026-09-16",
        "status": "SCOPED_FIXED_SPECTATOR_INSERTIONS_AND_UNIFORM_KNOWN_REFERENCE_RATE; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/insertion.md",
            "notes/cuts.md",
            "notes/bounds.md",
            "notes/matching.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_spectator_insertions": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_known_matching_and_uniform_rate": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The entire source-fixed H/Proca Gaussian metric kernels give complete crossed physical insertions with both optical cut normalizations and original M1 logs retained. Known finite shifts occupy the existing matching axes; remaining physical values stay unassigned. On the compact original domain the specified known-reference interference lies strictly between-3e-602 and0 at every nonforward angle, with endpoint-2dk/kappa. Gaussian-only disk and geometric bounds are not physical all-loop or Regge control. Original V/G/B/P8 remain OPEN.",
        "not_established": [
            "Values or bounds of the three physical finite matching coordinates",
            "Full detector-inclusive hard rate beyond the computed selected corrections",
            "Finite forward cross section or loop-square and all-order control",
            "Complex transfer analyticity, high-energy Regge and original quantum construction",
            "Original V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact conserved-stress projectors, both inverse-kernel signs, all crossed massive and massless cuts, source-fixed finite prescription and explicit written uniform bounds. Independent component contractions, spectral quadrature and complex-sheet calibrations supplement the general derivation. No unknown matching coordinate is chosen. Native/direct/ordinary/CLI use original SymPy; the adapter is FULL-only. Frozen ancestors, scoped P8(a) and original statuses are unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The spectator-insertion report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.305 fixed spectator insertion replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
