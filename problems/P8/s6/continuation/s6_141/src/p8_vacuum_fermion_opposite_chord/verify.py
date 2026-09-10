"""Read-only MS fermion opposite-chord certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_fermion_self_energy_chord import verify as parent

from . import angular, audit, calibration, catalog, domain, tail

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-fermion-opposite-chord.json"
PARENT_SHA = "f3355094e4a6a3a5a36692b56fb4ad263464e351726856fa8645419878f54fe6"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_fermion_opposite_chord/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen MS fermion opposite-chord report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_140_fully_rebuilt": PARENT_SHA,
        "same_GY14_nongravitational_reference_model": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
        "self_energy_subsets_not_complete_primitives_or_original_P8": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A MS fermion opposite-chord gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.141.SCALAR_GAUGE_OPPOSITE_CHORD_BOUNDS",
        "date": "2026-09-10",
        "status": "SCALAR_AND_GAUGE_OPPOSITE_CHORD_BOUNDS; NOT_COMPLETE_PRIMITIVE_ROWS_TWO_LOOP_MATCHING_HIGHER_LOOP_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/catalog.md",
            "notes/angular.md",
            "notes/tail.md",
            "notes/domain.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "joint_domain_and_exact_angular_integral": serialize(
            {"angular": payload(angular.data()), "domain": payload(domain.data())}
        ),
        "primitive_word_ownership": serialize(payload(catalog.data())),
        "soft_tail_and_radial_bound": serialize(payload(tail.data())),
        "actual_partial_primitive_frontier": serialize(
            {
                "actual": payload(calibration.data()),
                "primitive_frontier": audit.frontier(),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The twelve opposite-endpoint words in each scalar/gauge quartic fermion primitive have no divergent proper subgraph. Literal real-loop routing leaves the boson chord unshifted and gives a joint minimum-momentum Cauchy radius. Projecting the soft tail before both integrals permits an exact four-dimensional angular average and radial moment 2/3. The complete S4 subset has no lower-degree contribution to b2. The combined bound is below 1e-1400, relative to the tree below 1e-800. Together with S6.140, 36 of 60 words per sector are bounded; the 24 vertex words in each remain unbounded. Remaining rows, conversions, full two-loop error, V/G/B and original P8 remain open.",
        "not_established": [
            "Twenty-four vertex-correction words in each scalar/gauge quartic primitive",
            "Five other wholly unevaluated primitive fermion-sector rows",
            "Other counterterm, field and parameter-conversion contributions",
            "S6.138 finite outer-reference conversion to the shared MS boundary",
            "Complete two-loop amplitude/pole/error or higher-loop truncation control",
            "V contours, finite-gravity G, common-parent B or original P8 closure",
        ],
        "verification_boundary": "Exact word/cycle and routing checks, joint Dirac/Neumann domain identities, explicit angular and radial integrals, independent numerical checks and a mutation-checked partial frontier. Analytic regulator-limit and symmetry arguments are written proofs, not proof-assistant formalization or independent peer review. Native, ordinary, CLI and direct science use unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The MS fermion opposite-chord differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.141.SCALAR_GAUGE_OPPOSITE_CHORD_BOUNDS replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
