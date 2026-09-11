"""Read-only minimal-Einstein canonical-clock null-equation diagnostic."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_curved_dirac_stress import verify as parent

from . import audit, equations, map, nullstress

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-canonical-bounce-gap.json"
PARENT_SHA = "38ce293e50e209addb4bd4aa19f739be8a4dfc54ece5ac3f1890fc6bba449d60"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_canonical_bounce_gap/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen absolute curved free-stress report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_171_fully_rebuilt": PARENT_SHA,
        "same_actual_CD_free_in_out_Hadamard_states_and_EC_N0_reference": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A minimal canonical bounce-gap proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.172.MINIMAL_EINSTEIN_CANONICAL_CLOCK_NULL_EQUATION_OBSTRUCTION_AND_REQUIRED_LEADING_SOURCE",
        "date": "2026-09-11",
        "status": "EXCLUDED_NAMED_MINIMAL_EINSTEIN_CANONICAL_CLOCK_DIAGNOSTIC_WITH_SPECIFIED_FREE_STRESS_AND_SMALL_EXTRA_NULL_BUDGET; NOT_FULL_INTERACTING_PARENT_DHOST_ROW_V_G_B_OR_ORIGINAL_P8_EXCLUSION",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/model.md",
            "notes/null.md",
            "notes/map.md",
            "notes/budget.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "correlated_free_null_stress": serialize(payload(nullstress.data())),
        "minimal_parent_equation_and_map": serialize(
            {"equations": payload(equations.data()), "map": payload(map.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The full correlated free null stress is below1e595 uniformly, and below1e416 at the actual bounce after exact local cancellations. A minimally Einstein-coupled canonical clock instead leaves a positive null-equation residual at least9kappa minus that bounce allowance. On |t|<=1/2 the classical margin is at least121kappa/25. Additional classical positive-kinetic matter cannot repair it; neither can the unchanged clock-transparent map, whose first metric variation is identity even off shell. A repair at the bounce requires an additional null source <=-9kappa+the stated free allowance, not a perturbatively small correction to this diagnostic. No bound on the complete interacting remainder is derived, no original DHOST row is excluded, and V/G/B and original P8 remain open.",
        "not_established": [
            "An exclusion of the full GY14 interacting/gravitational completion, original CD action or any whole P8 ladder row",
            "A physical absolute bound on all additional quantum, nonminimal-gravity or higher-operator null stress",
            "Transfer of the free-state estimate to a nearby metric or arbitrary state; a controlled background or perturbation response",
            "Full matching, physical cutoff, vacuum contour/truncation or finite-gravity Regge/IR remainder",
        ],
        "verification_boundary": "Exact correlated local cancellations, physical-sign lapse/scale variations, positive-source algebra, clock-gate first variations and rational signed budgets support written proofs. Independent direct geometric curvature, finite variations, nonminimal-operator scope controls and complete-budget tests supplement the frozen free-state bounds. This is not a formalization or an interacting completion. Native, direct science, ordinary and CLI use original SymPy; only full regression uses its separately audited exact GCD adapter.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The minimal canonical bounce-gap report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.172.MINIMAL_EINSTEIN_CANONICAL_CLOCK_NULL_EQUATION_OBSTRUCTION_AND_REQUIRED_LEADING_SOURCE replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
