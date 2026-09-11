"""Read-only complete matched two-loop low-energy elastic-cut certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_two_loop_physical_source_map import verify as parent

from . import audit, bounds, calibration, cut, fermion, scalar

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-two-loop-elastic-cut.json"
PARENT_SHA = "e54116c62e54a0587bcb1e9dbb312284a905434a6db8276fdabef2842cd4dea5"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_two_loop_elastic_cut/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen low-energy cut parent report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_160_fully_rebuilt": PARENT_SHA,
        "same_GY14_canonical_physical_source_reference": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
        "finite_window_and_fixed_order_not_physical_remainder_or_global_V_G_B": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A complete low-energy cut proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.161.COMPLETE_MATCHED_TWO_LOOP_LOW_ENERGY_ELASTIC_CUT",
        "date": "2026-09-10",
        "status": "COMPLETE_FIXED_ORDER_TWO_LOOP_PHI_ELASTIC_CUT_ON_PHYSICAL_S_4_TO_6_AND_POSITIVE_FORMAL_IMPROVED_B2_MARGIN; NOT_PHYSICAL_TRUNCATION_GLOBAL_DISPERSION_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/ownership.md",
            "notes/scalar.md",
            "notes/fermion.md",
            "notes/cut.md",
            "notes/validation.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "full_physical_one_loop_amplitude_bounds": serialize(
            {
                "scalar": payload(scalar.data()),
                "fermion": payload(fermion.data()),
                "bounds": payload(bounds.data()),
            }
        ),
        "complete_two_loop_elastic_cut_and_ownership": serialize(payload(cut.data())),
        "actual_complete_two_loop_cut_and_improved_band": serialize(
            payload(calibration.data())
        ),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete first-loop canonical GY14 scattering amplitude on s in [4,6] is bounded using all full scalar MS bubbles/triangles/boxes, the full fermion local reference and convergent Dirac difference, and the sole first field correction. Resolving the scalar parameter boundary values before taking norms controls the light threshold. The only two-loop intermediate-state contribution on this window is the tree/one-loop two-Phi interference. Its integrated absolute allowance relative to tree b2 is below 1e-407. Subtracting both computed elastic-cut orders from the complete two-loop b2 polynomial preserves a positive formal margin with total relative allowance below 1e-6. Neither the second cut coefficient's sign nor the omitted physical remainder is assumed. Full dispersion, G, B and original P8 remain open.",
        "not_established": [
            "Physical higher-loop or operator-truncation errors for the amplitude, pole or cut",
            "A nonnegative exact spectral density inferred from the truncated interference polynomial",
            "Complete order-three gauge/four-particle cuts or high-energy contour/Regge control",
            "Ordinary-Psi composite normalization, finite-gravity G, common-parent B or original P8 closure",
        ],
        "verification_boundary": "Exact scalar parameter polynomials and boundary-weight derivatives, dimensional-reference-aware Dirac norm counting, full first-amplitude ownership, physical intermediate-state loop counting and rational cut/error composition. Independent tests compare log threshold integrals, boundary-subtracted physical triangles/boxes with analytic radial forms, direct Euclidean parameter integrals, explicit Dirac matrices and optical-theorem normalization controls. Source-pinned written analytic arguments and exact replay are not formalization or independent peer review. Native, direct science, ordinary and CLI retain unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete low-energy cut report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.161.COMPLETE_MATCHED_TWO_LOOP_LOW_ENERGY_ELASTIC_CUT replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
