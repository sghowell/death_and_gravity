"""Read-only finite forward phase and uniform real-interference certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_minimal_gravity_finite import verify as previous
from p8_vacuum_affine_one_newton_inclusive_assembly import verify as born_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-finite-forward-phase.json"
PARENT_SHA = "499d71347342e1d332367faf397757b617386dbf9faf4ebf17eb2be07b973780"
BORN_SHA = "6206f5d56959f01a8467e2c2b74ff43a2fad3cc75df0cee67e9bb8ce82f39564"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_finite_forward_phase/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen known finite forward-rate parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(born_input.REPORT) != BORN_SHA:
        raise ValueError("The frozen full positive Born normalization changed")
    born_input.validate_report(
        json.loads(born_input.REPORT.read_text()), born_input.build_report()
    )
    return {
        "S6_302_whole_known_finite_minimal_gravity_rebuilt": PARENT_SHA,
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
        raise ValueError("A finite forward-rate proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.303.COMPLETE_FINITE_FORWARD_COULOMB_PHASE_AND_UNIFORM_ALL_ANGLE_KNOWN_HARD_INTERFERENCE_BOUND",
        "date": "2026-09-15",
        "status": "SCOPED_ALL_ANGLE_KNOWN_REAL_INTERFERENCE_AND_FINITE_FORWARD_PHASE; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/phase.md",
            "notes/regularized.md",
            "notes/rate.md",
            "notes/matching.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_finite_forward_phase": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_regularized_forward_bound_and_real_rate": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The full finite Coulomb phase and real classical/logarithmic forward terms of the unchanged known minimal-gravity reference are explicit. Exact pole cancellation and uniform derivative bounds give a known real hard-reference one-loop interference below1e12/kappa, hence1e-788 at originalparameters, over every physical nonforward angle in the compact energy range. The complex logarithmic phase remains. Unknown matching, complete radiation and all-loop/Regge requirements and original V/G/B/P8 remain OPEN.",
        "not_established": [
            "Values or bounds of the three physical finite matching coordinates",
            "Gravity-Born radiation and full detector-inclusive rate",
            "Finite forward cross section or loop-square and all-order control",
            "Complex transfer analyticity, high-energy Regge and original quantum construction",
            "Original V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact forward coefficients, full massive-block cancellation, weighted rational derivative norms and Born inequalities accompany written endpoint and uniform-rate proofs. Independent high-precision complete-master calculations calibrate both endpoint coefficients and all-angle rates, not rigorous quadrature errors. The imaginary phase is retained in the complex amplitude. Native/direct/ordinary/CLI use original SymPy; the adapter is FULL-only. Frozen ancestors and original statuses are unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The finite forward-rate report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.303 finite forward-rate reference replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
