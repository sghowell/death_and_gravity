"""Read-only replay for a separate regular classical affine/vector parent."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_target_clock_null import verify as parent

from . import affine, audit, chart, dynamics, source

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-analytic-affine-parent.json"
PARENT_SHA = "3740dfeae08879dda0039c10fd94c09018ef964d1ca12690624089927a852ba5"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_analytic_affine_parent/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen leading target-clock report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_173_fully_rebuilt": PARENT_SHA,
        "same_complete_analytic_target_and_original_M1": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A regular classical-parent proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.174.REGULAR_CLASSICAL_SOURCE_ALIGNED_AFFINE_PARENT_WITH_EXACT_CD_CLOCK_AND_VACUUM_QUADRATIC_BLOCKS",
        "date": "2026-09-11",
        "status": "SEPARATE_REGULAR_CLASSICAL_AFFINE_VECTOR_PARENT_EXACT_TARGET_REDUCTION_AND_LOCAL_CLOCK_CONSTRAINTS; NOT_QUANTUM_MATCHING_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/chart.md",
            "notes/affine.md",
            "notes/source.md",
            "notes/constraints.md",
            "notes/modes.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "regular_classical_parent_and_full_algebraic_reduction": serialize(
            {
                "chart": payload(chart.data()),
                "source": payload(source.data()),
                "affine": payload(affine.data()),
            }
        ),
        "actual_clock_vacuum_and_constraint_blocks": serialize(
            {
                "target": payload(chart.target_jets()),
                "dynamics": payload(dynamics.data()),
                "clock": payload(dynamics.clock_lapse()),
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The separately named CD-REG-AFFINE-ISO classical EFT retains the complete analytic S6.109 target and original M1. A regular covariant metric chart, the vacuum-regular lower dictionary, and a new source-centered isotropic trace mass give an exact full-quotient reduction to that target plus a retained Proca field with nonlinear source. The source and its first variation vanish on the entire CD clock, giving an exact classical clock and the unchanged light quadratic block plus positive Proca. The joint lapse/vector auxiliary Jacobian is nonsingular at every finite clock point, with local seven-mode neighborhoods; the vacuum has the stated healthy seven-mode quadratic spectrum. This constructs classical common-domain structure, not quantum target matching, a physical cutoff, a controlled interacting state/background response, or V/G/B and original P8 closure.",
        "not_established": [
            "An interacting quantum parent, transfer of GY14 matching or free-sector state/stress bounds",
            "A physical cutoff, all omitted-order bounds, or uniform nonlinear stability/causality on the connected domain",
            "Controlled vector integration, quantum-corrected background or perturbation response",
            "Vacuum contour/truncation, finite-gravity IR/Regge remainder, V/G/B or original P8 closure",
        ],
        "verification_boundary": "Full source-pinned generic-R covariant/ADM identities, all60 quotient Euler lifts, source centering, explicit analytic inverse/limits and exact clock lapse polynomial support written proofs. Independent tensor and finite metric variations, dense quotient/constraint solves, exact spectral limits and numerical full-function controls supplement them. Local constraint counting is not a formalized theorem or a nonlinear stability claim. Native, direct science, ordinary and CLI use original SymPy; full regression alone uses its separately audited exact GCD adapter.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The regular analytic affine-parent report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.174.REGULAR_CLASSICAL_SOURCE_ALIGNED_AFFINE_PARENT_WITH_EXACT_CD_CLOCK_AND_VACUUM_QUADRATIC_BLOCKS replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
