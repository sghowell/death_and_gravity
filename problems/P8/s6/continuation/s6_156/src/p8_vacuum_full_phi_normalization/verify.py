"""Read-only complete GY14 two-loop Phi normalization and local pole certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_finite_field_covariance import verify as parent

from . import audit, bounds, calibration, normalization, ownership, pole

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-full-phi-normalization.json"
PARENT_SHA = "76f6fac15d1f81af93b428f2d074ea12d12f5bbfcf08a6414a9dd2035822ed82"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_full_phi_normalization/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen complete Phi normalization parent report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_155_fully_rebuilt": PARENT_SHA,
        "same_GY14_nongravitational_reference_model": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
        "complete_Phi_field_and_local_pole_advance_not_full_amplitude_vacuum_or_V_G_B": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A complete Phi normalization proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.156.COMPLETE_GY14_TWO_LOOP_PHI_NORMALIZATION_AND_POLE",
        "date": "2026-09-10",
        "status": "COMPLETE_GY14_FIXED_ORDER_TWO_LOOP_PHI_NORMALIZATION_AND_CANONICAL_UNIT_DISC_POLE; NOT_FULL_FOUR_POINT_VACUUM_SOURCE_TRUNCATION_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/ownership.md",
            "notes/normalization.md",
            "notes/matching.md",
            "notes/pole.md",
            "notes/bounds.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "regulator_complete_field_and_scalar_parameter_map": serialize(
            payload(normalization.data())
        ),
        "complete_quadratic_ownership_and_OS_extension": serialize(
            {
                "ownership": payload(ownership.data()),
                "pole": payload(pole.data()),
                "bounds": payload(bounds.data()),
            }
        ),
        "actual_complete_Phi_two_loop_enclosure": serialize(
            payload(calibration.data())
        ),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "All 32 scalar refinements and all three fermionic quadratic families, with their assigned proper counterterms and physical inner OS grouping, give a complete hybrid two-loop inverse slope. The full regulated first-parameter shift cancels the epsilon-times-MS-field-pole product, yielding t_MS=t_H-k0 r0. Re-expressing the first interaction poles also fixes the second scalar G,L,M parameter map, including the induced heavy-mass shift. The complete MS second Phi normalization is bounded below 1e-18. Extending the direct primitive Cauchy estimate to the unit disc and retaining both first-sector parameter variations bounds the full canonical fixed-order inverse by |Pi_R|<1e-18 |s-1|^2. It therefore has only the mass-one unit-residue light pole on that disc. This is not the matched four-point/vacuum/source assembly, a global or all-orders spectrum, physical truncation or V/G/B closure.",
        "not_established": [
            "The complete matched four-point amplitude and vacuum/source reference assembly",
            "External-fermion second-order canonical dictionaries beyond the Phi calculation",
            "Global spectral/quantum-potential results or physical higher-loop truncation",
            "V contours/cuts, finite-gravity G, common-parent B or original P8 closure",
        ],
        "verification_boundary": "Complete explicit quadratic ownership, common-regulator field and bare-coordinate identities, inherited integrated first-sheet bounds and exact rational enclosures. Tests use independent regulated coupling substitutions, finite Laurent extraction, Cauchy-tail examples and pole/normalization mutation controls. Written analytic arguments and exact replay are not formalization or independent peer review. Native, ordinary, CLI and direct science retain unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete Phi normalization report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.156.COMPLETE_GY14_TWO_LOOP_PHI_NORMALIZATION_AND_POLE replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
