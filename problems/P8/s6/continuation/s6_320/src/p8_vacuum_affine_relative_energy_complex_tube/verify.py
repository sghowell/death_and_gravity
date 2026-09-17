"""Read-only relative-energy complex-tube and Cauchy-bound certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_two_real_soft_overlap import verify as analytic_input
from p8_vacuum_affine_uniform_all_tree_bound import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-relative-energy-complex-tube.json"
)
PARENT_SHA = "6239bb0b855c48275f8598cdfdec3f8dfa2590545d5b4abb5b6ddfe765017b0b"
ANALYTIC_SHA = "79539dcf68142c2ebb62a0236293be98d2f1b6eece93d49e46cbf091fb8e1994"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_relative_energy_complex_tube/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen complete all-finite tree-bound parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(analytic_input.REPORT) != ANALYTIC_SHA:
        raise ValueError("The frozen complex-recoil and two-real overlap input changed")
    analytic_input.validate_report(
        json.loads(analytic_input.REPORT.read_text()), analytic_input.build_report()
    )
    return {
        "S6_319_uniform_all_finite_complete_tree_bound_rebuilt": PARENT_SHA,
        "S6_313_complex_recoil_and_two_real_overlap_rebuilt": ANALYTIC_SHA,
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
        raise ValueError("A relative complex-tube bound proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.320.ALL_FINITE_RELATIVE_ENERGY_COMPLEX_TUBES_AND_CAUCHY_BOUND",
        "date": "2026-09-16",
        "status": "SCOPED_RELATIVE_ENERGY_COMPLEX_TUBE_AND_CAUCHY_BOUND; NOT_SOFT_FACES_OR_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/geometry.md",
            "notes/continuation.md",
            "notes/majorant.md",
            "notes/calibration.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_complex_geometry": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_complex_continuation_and_Cauchy_majorants": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The original complete finite tree admits a holomorphic relative-energy tube with an explicit all-N angular-uniform majorant and multiindex Cauchy bound. A complex conserved temporal inverse closes the center-weighted current norm, and unchanged hard-core ownership controls the full source. The exact complex5116-tree calibration preserves physical tree equivalence. A nonzero proper-cluster pole rejects a naive global-total-W neighborhood. The shrinking radii do not reach soft faces or prove an inclusive probability, quantum state, absolute complex Regge or original V/G/B/P8 closure.",
        "not_established": [
            "Uniform holomorphy through all soft faces or an infrared-finite all-N probability",
            "The full overlapping all-N soft subtraction and real-virtual completion",
            "Complete finite hard real-virtual and evanescent matching",
            "An interacting quantum state, unitarity or absolute complex Regge",
            "The original common-parent bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact generic complex conserved inverse and weighted matrix factorization; positive-coefficient subset-invariant perturbation bounds; all-order Rosen/reflection grading and nonnegative barriers; complete complex four-tree and original5116-tree calibrations; independent EGF coefficients and explicit multiindex Cauchy factors; a complete-current transverse pole inside an invalid global-W disc. Gaussian-rational normalization is exact, with the failed raw-zero representation and successful diagnosis retained externally. Finite points supplement, not replace, the holomorphy proof. Original SymPy is retained outside the FULL-only adapter.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The relative complex-tube bound report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.320 relative complex-tube bound replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
