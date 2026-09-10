"""Read-only scalar and gauge fermion vacuum primitive certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_fermion_insertion_ms import verify as parent

from . import audit, calibration, derivatives, forest, massless, scalar

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-fermion-vacuum.json"
PARENT_SHA = "b7740228c719dede347afa67d3b5853c677c984ba4a15bfa4209836e38c8408d"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_fermion_vacuum/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen vacuum primitive parent report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_148_fully_rebuilt": PARENT_SHA,
        "same_GY14_nongravitational_reference_model": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
        "only_two_vacuum_primitive_rows_advance": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A vacuum primitive proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.149.PAIRED_FERMION_VACUUM_PRIMITIVES",
        "date": "2026-09-10",
        "status": "BOTH_VACUUM_PRIMITIVES_BOUNDED_WITH_ASSIGNED_PROPER_FORESTS; ALL_NINE_PRIMITIVE_ROWS_HAVE_SCOPED_BOUNDS; NOT_FULL_VACUUM_MATCHING_CANONICAL_POLE_TRUNCATION_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/massless.md",
            "notes/scalar.md",
            "notes/forests.md",
            "notes/remainder.md",
            "notes/reference.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "exact_massless_vacuum_and_derivative_anchors": serialize(
            {
                "massless": payload(massless.data()),
                "derivatives": payload(derivatives.data()),
            }
        ),
        "actual_scalar_vacuum_and_assigned_forests": serialize(
            {"scalar": payload(scalar.data()), "forest": payload(forest.data())}
        ),
        "actual_vacuum_primitive_enclosure": serialize(payload(calibration.data())),
        "partial_primitive_frontier": serialize(audit.frontier()),
        "controls": serialize(audit.controls()),
        "verdict": "The two vacuum primitive rows retain all assigned proper forests and their finite outer MS references. Dimension-symbolic vacuum tensors and proper fermion counterterms give massless-exchange finite coefficients -19 NYm^4/Q^2 and +18 NaC_Fm^4/Q^2. The actual scalar row is paired with its complete physical mass/residue subgraph and evaluated through three beta anchors and a bounded finite remainder. Gauge exchange includes all fourteen Dirac flavors. The combined absolute vacuum contribution is below 1e595 and below 1e-203 relative to the one-loop fermion vacuum reference. All nine primitive rows now have bounds in their stated scopes, not a complete matching or P8 verdict. Other vacuum/source contributions, full canonical matching, truncation and V/G/B remain open.",
        "not_established": [
            "Complete whole-model vacuum-reference counterterm including other scalar and source terms",
            "Other parameter/counterterm matching insertions and complete canonical two-loop amplitude",
            "Complete matched pole/error, higher-loop truncation or V contour/cut control",
            "Finite-gravity G, common-parent B or original P8 closure",
        ],
        "verification_boundary": "Exact dimension-symbolic vacuum traces, full assigned proper forests, fixed-scale derivative checks, Gamma duplication/recurrences, three beta anchors and an endpoint-controlled finite remainder. Independent checks use explicit Dirac matrices, nonzero-regulator integrals and finite-part extraction. Written analytic arguments and exact replay are not formalization or independent peer review. Native, ordinary, CLI and direct science retain unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The vacuum primitive report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.149.PAIRED_FERMION_VACUUM_PRIMITIVES replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
