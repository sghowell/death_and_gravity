"""Read-only certificate for the two quadratic primitive finite MS masses."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_fermion_ms_slopes import verify as parent

from . import anchors, audit, calibration, correction, forest, outer

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-fermion-ms-mass.json"
PARENT_SHA = "8df7827728af867a091e4532574132d530f0e3815f4473cb9404e918c27ae202"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_fermion_ms_mass/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen finite MS mass parent report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_146_fully_rebuilt": PARENT_SHA,
        "same_GY14_nongravitational_reference_model": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
        "only_two_primitive_finite_MS_masses_advance": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A finite MS mass proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.147.FINITE_MS_QUADRATIC_PRIMITIVE_MASSES",
        "date": "2026-09-10",
        "status": "BOTH_QUADRATIC_PRIMITIVE_NONLOCAL_AND_FINITE_MS_MASS_AND_SLOPE_BOUNDS; NOT_OTHER_FORESTS_FULL_MATCHING_CANONICAL_POLE_HIGHER_LOOPS_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/anchors.md",
            "notes/forests.md",
            "notes/mass_integral.md",
            "notes/remainder.md",
            "notes/pole.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "exact_massless_MS_anchors": serialize(payload(anchors.data())),
        "full_regulated_scalar_mass_transfer": serialize(
            {
                "forest": payload(forest.data()),
                "outer": payload(outer.data()),
                "correction": payload(correction.data()),
            }
        ),
        "actual_primitive_mass_enclosure": serialize(payload(calibration.data())),
        "partial_primitive_frontier": serialize(audit.frontier()),
        "controls": serialize(audit.controls()),
        "verdict": "The complete dimensional quadratic primitive mass tensors and proper fermion counterterms give finite massless-exchange references -56 NY^2 m^2/Q^2 and +40 NYaC_F m^2/Q^2. The actual scalar mass is restored by integrating the complete regulated mixed-quartic outer reference, including its whole-cycle quartic MS counterterm. Its leading finite contribution is NY^2/Q^2 [78+32 log(m^2)], with an explicit finite-ratio remainder. Previous slope and nonlocal bounds give a combined on-shell absolute mass-reference bound below 1e-10. Other forests, the older W1F0 local references, two vacuum rows, full matching and original P8 remain open.",
        "not_established": [
            "Finite MS local references of the older scalar_Phi2_W1_F0 row",
            "Two vacuum primitive rows and other counterterm/matching insertions",
            "Complete matched/canonical two-loop amplitude, pole or error",
            "Higher-loop truncation control, V contours, finite-gravity G, common-parent B or original P8 closure",
        ],
        "verification_boundary": "Exact dimensional mass tensors, independent vacuum differentiation, full proper-counterterm Laurent expansion and regulated scalar-mass integration. Numerical tests independently extract Laurent coefficients, integrate general-mass triangles and finite remainders, and check explicit Dirac matrices. Analytic arguments are written proofs, not formalization or independent peer review. Native, ordinary, CLI and direct science retain unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The finite MS mass report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.147.FINITE_MS_QUADRATIC_PRIMITIVE_MASSES replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
