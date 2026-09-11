"""Read-only complete stationary-H source and fixed-order reference certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_full_vacuum_reference import verify as parent

from . import audit, bounds, calibration, differentiation, ownership, source

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-full-heavy-source.json"
PARENT_SHA = "5867f46444e1227ea22cce2cfc70ac5fc4c8e28893e90ad5c2ae1e9632b22430"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_full_heavy_source/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen complete H-source parent report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_158_fully_rebuilt": PARENT_SHA,
        "same_GY14_nongravitational_reference_model": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
        "named_fixed_order_reference_assembly_not_other_dictionaries_truncation_or_V_G_B": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A complete H-source proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.159.COMPLETE_SECOND_H_SOURCE_AND_FIXED_ORDER_REFERENCE_ASSEMBLY",
        "date": "2026-09-10",
        "status": "COMPLETE_MATCHED_GY14_TWO_LOOP_PHI_POLE_B2_VACUUM_AND_STATIONARY_H_SOURCE_REFERENCES; NOT_OTHER_COORDINATE_DICTIONARIES_GLOBAL_POTENTIAL_PHYSICAL_TRUNCATION_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/ownership.md",
            "notes/source.md",
            "notes/matching.md",
            "notes/derivative.md",
            "notes/bounds.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_stationary_source_ownership_and_matching": serialize(
            {"ownership": payload(ownership.data()), "source": payload(source.data())}
        ),
        "mass_differentiated_sunset_and_source_bounds": serialize(
            {
                "differentiation": payload(differentiation.data()),
                "bounds": payload(bounds.data()),
            }
        ),
        "actual_complete_second_H_source_enclosure": serialize(
            payload(calibration.data())
        ),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete second stationary-H source contains the scalar and fermionic OS covariance tadpoles, the proper cubic MS counterterm tadpole and the full first-Phi-coordinate re-expression of the regulated first source. Differentiating the vacuum functional with fixed reference counterterms gives the same scalar insertion; varying those counterterms with H would be incorrect. Six differentiated sunset sectors and an explicit small-r split bound the remaining scalar terms. The positive-epsilon first field coefficient times the source pole is retained. The complete source bound is below 1e189, its normalized H reference below 1e-8 and its induced Phi mass reference below 1e-12. Together with the unchanged complete Phi pole, b2 and vacuum reference, this completes the named fixed-order reference assembly. Other coordinate dictionaries, global potential, physical truncation and V/G/B remain open.",
        "not_established": [
            "External-fermion or derivative-coordinate/source dictionaries beyond the named Phi/vacuum reference",
            "A global effective potential, nonperturbative vacuum or semiclassical backreaction theorem",
            "A physical finite-EFT higher-order truncation bound or higher-order regulated source-square reference",
            "V contours/cuts, finite-gravity G, common-parent B or original P8 closure",
        ],
        "verification_boundary": "Complete one-point counterterm ownership, fixed-reference vacuum differentiation, full regulated bare-source re-expression, six mass-differentiated sectors and exact rational power-cap bounds. Independent tests compare OS covariance and fixed-counterterm vacuum derivatives, reject moving reference conditions, extract finite source products on complex regulator circles and check the two-light-line sunset mass derivative and endpoint bounds. Written analytic arguments and exact replay are not formalization or independent peer review. Native, ordinary, CLI and direct science retain unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The complete H-source report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.159.COMPLETE_SECOND_H_SOURCE_AND_FIXED_ORDER_REFERENCE_ASSEMBLY replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
