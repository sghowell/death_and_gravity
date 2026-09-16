"""Read-only full first-Newton inclusive-assembly certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_heavy_parent_one_loop import verify as matter_input
from p8_vacuum_affine_physical_virtual_soft_pairing import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-one-newton-inclusive-assembly.json"
)
PARENT_SHA = "40c04cadfa1a1c542c16eed48f6df12bcfe0605ad69b39dbb70fcefaefe6b3ed"
MATTER_SHA = "a09852126d79ca223af06555e10283328cb144a14945c127abee73b7ce8b9c05"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_one_newton_inclusive_assembly/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen physical soft-pairing parent report changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(matter_input.REPORT) != MATTER_SHA:
        raise ValueError("The frozen full matter loop changed")
    matter_input.validate_report(
        json.loads(matter_input.REPORT.read_text()), matter_input.build_report()
    )
    return {
        "S6_296_physical_soft_pairing_entire_ancestry_rebuilt": PARENT_SHA,
        "S6_239_complete_same_OS4_matter_loop_rebuilt": MATTER_SHA,
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
        raise ValueError("A complete first-Newton assembly proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.297.COMPLETE_FORMAL_ONE_NEWTON_INCLUSIVE_ASSEMBLY_FINITE_KAPPA_MATTER_GRAPH_EQUALITY_AND_FORWARD_CROSSOVER",
        "date": "2026-09-15",
        "status": "SCOPED_COMPLETE_FORMAL_ONE_NEWTON_INCLUSIVE_ASSEMBLY; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/matter.md",
            "notes/forward.md",
            "notes/inclusive.md",
            "notes/matching.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_actual_finite_kappa_source_and_complete_matter_loop": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_positive_Born_forward_boundary_and_inclusive_assembly": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The actual finite-kappa source has the same complete contributing matter vertices and first-loop graphs as S239 in its existing OS4 prescription. The missing gravity-Born times matter-loop rate interference is uniformly below10^-199 after positive full Born normalization. Combined with S296, the entire formal first-Newton one-loop inclusive numerator is assembled, with the unknown gravitational hard matching explicit and a reference error below2*10^-199 on the compact detector domain. An exact forward crossover rules out an all-angle small-Newton inference. No full forward, all-energy Regge, omitted-loop or original V/G/B/P8 closure is claimed.",
        "not_established": [
            "Independent finite gravitational matching or a UV-completion verdict",
            "Higher-Newton or omitted-loop error and full forward cross section",
            "A uniform forward Newton expansion or exchanged infrared limits",
            "Crossing-analytic positive dispersion measure or bounded Regge remainder",
            "Original quantum state/domain/measure/bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Whole finite-kappa vertices, exact graph counting, the fully rebuilt earlier matter loop, and positive Born/error identities establish this formal perturbative assembly. It neither commutes quantization with a limit nor declares the inclusive rate an analytic amplitude. The full Born normalization is explicit and does not manufacture missing higher Newton orders. The finite-contour discussion is a conditional boundary statement, not a verified model-specific Regge estimate. Native/direct/ordinary/CLI retain original SymPy; the exact-GCD adapter remains FULL-only. Frozen scientific source/raw evidence and original closure statuses are unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete first-Newton inclusive report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.297 full first-Newton inclusive assembly replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
