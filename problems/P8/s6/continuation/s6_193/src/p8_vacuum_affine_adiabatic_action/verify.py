"""Read-only actual fixed covariant shear matching response and finite-amplitude remainder."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_current_response import verify as parent
from p8_vacuum_affine_proca_gaussian import verify as gaussian

from . import action, angular, audit, curvature, matching

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-adiabatic-action.json"
PARENT_SHA = "3e7a44038928cdcb613d07545ac27418fc0dc21d819d19e9407692f57ff2dba7"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_adiabatic_action/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    for previous, expected in ((parent, PARENT_SHA), (gaussian, GAUSSIAN_SHA)):
        if sha(previous.REPORT) != expected:
            raise ValueError(
                "The source-pinned actual affine Gaussian state chain changed"
            )
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_192_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError(
            "An actual fixed covariant shear matching response gate failed"
        )
    return {
        "schema": 1,
        "claim": "P8-S6.193.ACTUAL_HOMOGENEOUS_SHEAR_FULL_NONCOMMUTING_ADIABATIC_CURRENT_ACTION_AND_FIXED_COVARIANT_SUBTRACTION_MATCHING",
        "date": "2026-09-11",
        "status": "ACTUAL_HOMOGENEOUS_SHEAR_FIXED_COVARIANT_SUBTRACTION_AND_FINITE_CONTACT_MATCHING; NOT_FULL_SPATIAL_MIXED_RESPONSE_INVERSE_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/action.md",
            "notes/variation.md",
            "notes/angular.md",
            "notes/dimension.md",
            "notes/matching.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "full_local_action_and_covariant_matching": serialize(
            {
                "action": payload(action.data()),
                "angular": payload(angular.data()),
                "matching": payload(matching.data()),
            }
        ),
        "dimension_limit_and_fixed_prescription": serialize(payload(curvature.data())),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The actual fourth-order homogeneous-shear mode current is the compact variation of the full local finite marker action, as verified for every noncommuting matrix word and independent frequency trace. Complete oriented angular contractions and radial continuation in general dimension give the covariant curvature expression modulo explicit compact total derivatives. Keeping the original four-dimensional scalar pole coefficients and all evanescent polarization and curvature factors reproduces exactly the S176 finite local density, including 5/2,5/3 and the -4a4s contact. The original fixed-prescription homogeneous current is therefore the convergent S192 comparison plus the variation of that unchanged local action. No state, parent or counterterm is selected anew. A numerical nonlinear local-current norm and full spatial/mixed response, feedback inverse, interacting finite-coupling background, physical cutoff and original V/G/B and P8 closure are not established.",
        "not_established": [
            "A numerical nonlinear local-current norm or complete spatial/mixed physical response",
            "A full feedback inverse, finite-coupling interacting quantum background or its stability",
            "Full parent state/measure/loops, physical heavy/cutoff/threshold and omitted-order matching",
            "Complete vacuum cuts/contour/truncation, finite-gravity IR/Regge or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact full noncommuting word/current variations, sphere pairings, general-d shear curvature coefficients, compact boundary identities and finite Laurent coefficients. Independent Cartesian sphere integration, non-diagonal Christoffel/Riemann and literal metric/frequency Euler fixtures supplement written proofs. Native, direct science, ordinary and CLI use original SymPy; only full regression uses the audited exact GCD adapter. Analytic and variational arguments are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The actual fixed covariant shear matching report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.193.ACTUAL_HOMOGENEOUS_SHEAR_FULL_NONCOMMUTING_ADIABATIC_CURRENT_ACTION_AND_FIXED_COVARIANT_SUBTRACTION_MATCHING replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
