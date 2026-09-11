"""Read-only leading actual CD clock sources and matching-domain audit."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_canonical_bounce_gap import verify as parent

from . import audit, homogeneous, structural, target

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-target-clock-null.json"
PARENT_SHA = "06a5393e97b58301c4b503dd20ff3758158d91cad1dd59f8926750ee274fcaa5"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_target_clock_null/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen minimal canonical bounce-gap report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_172_fully_rebuilt": PARENT_SHA,
        "same_actual_CD_target_clock_and_original_rolling_matter": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("An actual target-clock decomposition proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.173.ACTUAL_CD_LEADING_CLOCK_EULER_SOURCES_AND_CONSTRAINT_MATCHING_DOMAIN_GAP",
        "date": "2026-09-11",
        "status": "EXACT_LEADING_ACTUAL_CD_TARGET_BACKGROUND_SOURCE_DECOMPOSITION_WITH_CONSTRAINT_AND_DOMAIN_CONTROLS; NOT_PROPAGATING_PARENT_QUANTUM_BOUNCE_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/geometry.md",
            "notes/target.md",
            "notes/matching.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_covariant_variations_and_structural_controls": serialize(
            {
                "homogeneous": payload(homogeneous.data()),
                "structural": payload(structural.data()),
            }
        ),
        "actual_leading_clock_source_decomposition": serialize(payload(target.data())),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged full analytic CD_matter target retains the exact first clock jets. Varying the covariant action before fixing the lapse gives curvature and A3 null sources+24kappa and-12kappa at the bounce, while the scalar correction relative to a canonical clock gives-2101kappa/100. The original M1 source+kappa/100 is retained; both full metric equations cancel exactly at every time. A4,A5 have zero first background variation but are required for the clock velocity-Hessian degeneracy. The unit Fourier-jet vacuum matching class remains separated from the canonical clock by a normalized first-jet gap above1/2. These identify leading common-parent requirements; they do not construct the parent, control quantum-corrected background response, or close V/G/B and original P8.",
        "not_established": [
            "A propagating healthy UV parent reproducing the full leading clock action and constraint neighborhood",
            "Quantum target matching, interacting state/cutoff, omitted-order bounds or state-aware background/perturbation response",
            "An extension of the vacuum Fourier-unit action estimate to the large-gradient unit clock",
            "Vacuum contour/truncation, finite-gravity Regge/IR remainder, V/G/B or original P8 closure",
        ],
        "verification_boundary": "Exact full-target first clock jets, complete homogeneous covariant first metric variations, both entire-time equations, original matter conservation and trace/lapse Hessian controls support written proofs. Independent direct tensor contractions, functional metric variations, matter Routhian/sign and matching-domain controls supplement them. These do not formalize continuum theorems or certify an interacting completion. Native, direct science, ordinary and CLI use original SymPy; only full regression uses its separately audited exact GCD adapter.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The actual target-clock null report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.173.ACTUAL_CD_LEADING_CLOCK_EULER_SOURCES_AND_CONSTRAINT_MATCHING_DOMAIN_GAP replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
