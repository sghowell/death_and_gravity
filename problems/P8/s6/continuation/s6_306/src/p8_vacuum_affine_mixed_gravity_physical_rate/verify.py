"""Read-only physical mixed-gravity hard-reference rate certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_one_newton_inclusive_assembly import verify as born_input
from p8_vacuum_affine_spectator_gravity_insertion import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-mixed-gravity-physical-rate.json"
PARENT_SHA = "22b7a327ff4d4a2c39fc6ef4fe7435e84549def254317347e7f8abc01a362011"
BORN_SHA = "6206f5d56959f01a8467e2c2b74ff43a2fad3cc75df0cee67e9bb8ce82f39564"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_mixed_gravity_physical_rate/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen complete spectator insertion parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(born_input.REPORT) != BORN_SHA:
        raise ValueError("The frozen full positive Born normalization changed")
    born_input.validate_report(
        json.loads(born_input.REPORT.read_text()), born_input.build_report()
    )
    return {
        "S6_305_whole_fixed_spectator_insertions_and_known_rate_rebuilt": PARENT_SHA,
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
        raise ValueError("A mixed-gravity-rate proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.306.PHYSICAL_MIXED_GRAVITY_KNOWN_RATE_AND_EXACT_WHOLE_BOX_COLLAPSE",
        "date": "2026-09-16",
        "status": "SCOPED_COMPLETE_PHYSICAL_MIXED_GRAVITY_KNOWN_RATE; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/masters.md",
            "notes/contour.md",
            "notes/bounds.md",
            "notes/assembly.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_physical_masters": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_finite_assembly_and_uniform_rate": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete known S293/S294/S290 mixed-sector graph sum now has an explicit physical-sheet hard-reference bound. An exact projective identity collapses all six mixed boxes. Full-D quartic finite terms, both massive endpoint triangles, OS terms and the H-metric bubble remain, without C-endpoint double counting. At original parameters the linear interference normalized by the full positive Born is below4e-204 uniformly at all nonforward angles and tends to0 at either angular endpoint. Unknown physical matching and original V/G/B/P8 remain OPEN.",
        "not_established": [
            "Values or bounds of finite physical curvature and hard matching",
            "Uncomputed loop squares and complete all-order detector-inclusive rate",
            "Finite forward cross section or exchanged regulator and resolution limits",
            "Complex transfer and high-energy Regge estimates or a UV quantum construction",
            "Original V/G/B/P8 closure",
        ],
        "verification_boundary": "Whole inherited graph identity, exact projective Jacobian, explicit nonpinching physical contour, full radial endpoint moments and positive arithmetic bounds. Independent original-parameter quadrature and the primary Box16 formula calibrate normal-sheet finite terms; numerical samples are not the uniform proof. Native/direct/ordinary/CLI use original SymPy; the adapter is FULL-only. Frozen ancestors, scoped P8(a), matching coordinates and original statuses are unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The mixed-gravity-rate report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.306 physical mixed-gravity known-rate replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
