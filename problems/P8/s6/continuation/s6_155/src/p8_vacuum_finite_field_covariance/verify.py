"""Read-only regulated first finite Phi field covariance certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_scalar_ms_slopes import verify as parent

from . import audit, bounds, calibration, coefficients, covariance, kernel

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-finite-field-covariance.json"
PARENT_SHA = "67cd5c13683398935e985722d2018966d3bfcd25ca0f7657d949fdf619c08a0e"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_finite_field_covariance/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen finite Phi field parent report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_154_fully_rebuilt": PARENT_SHA,
        "same_GY14_nongravitational_reference_model": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
        "only_first_finite_Phi_normalization_direction_advances": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A finite Phi field proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.155.REGULATED_FIRST_FINITE_PHI_FIELD_COVARIANCE",
        "date": "2026-09-10",
        "status": "EXACT_REGULATED_FIRST_FINITE_PHI_FIELD_DIRECTION_AND_BOUNDED_ORDER_TWO_CROSS_TERMS; NOT_NEW_SECOND_SLOPE_FULL_GY14_MATCHING_TRUNCATION_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/kernel.md",
            "notes/covariance.md",
            "notes/coefficients.md",
            "notes/amplitude.md",
            "notes/bounds.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "full_regulated_first_slope_and_MS_coefficient_map": serialize(
            {
                "kernel": payload(kernel.data()),
                "coefficients": payload(coefficients.data()),
            }
        ),
        "field_covariance_and_first_slope_enclosures": serialize(
            {"covariance": payload(covariance.data()), "bounds": payload(bounds.data())}
        ),
        "actual_first_field_order_two_and_coordinate_bounds": serialize(
            payload(calibration.data())
        ),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The first finite Phi normalization direction is organized through formal order two at the common regulator. Both scalar and fermion slope epsilon coefficients are determined; the latter follows a full dimensional trace reduction and parameter integration by parts. Regulated vertices and internal lines combine to a single external field factor, including assigned counterterms. The minimal-pole coefficient map retains its finite epsilon-times-pole commutator and the fundamental G and y square cross terms. That coordinate term is not added again to the fully renormalized amplitude. The second-order amplitude terms determined solely by the first normalization are bounded below 1e-212 relative to tree. The genuinely new second slope remains symbolic, and other matching directions, full GY14 assembly, physical truncation and V/G/B remain open.",
        "not_established": [
            "The genuinely new complete two-loop scalar field normalization coefficient",
            "Other finite parameter directions, their cross terms and full GY14 field/coupling assembly",
            "Complete vacuum/source and canonical amplitude/pole assembly or physical higher-loop error",
            "V contours/cuts, finite-gravity G, common-parent B or original P8 closure",
        ],
        "verification_boundary": "Exact full-dimensional slope/IBP algebra, common-regulator vertex/line covariance, Laurent coefficient maps and inherited complete one-loop enclosures. Tests independently evaluate nonzero-regulator trace derivatives, complex finite-part coefficients, literal field/propagator transformations and nonlinear G/y expansions, with input and scope mutations. Written proofs and exact replay are not formalization or independent peer review. Native, ordinary, CLI and direct science retain unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The finite Phi field report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.155.REGULATED_FIRST_FINITE_PHI_FIELD_COVARIANCE replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
