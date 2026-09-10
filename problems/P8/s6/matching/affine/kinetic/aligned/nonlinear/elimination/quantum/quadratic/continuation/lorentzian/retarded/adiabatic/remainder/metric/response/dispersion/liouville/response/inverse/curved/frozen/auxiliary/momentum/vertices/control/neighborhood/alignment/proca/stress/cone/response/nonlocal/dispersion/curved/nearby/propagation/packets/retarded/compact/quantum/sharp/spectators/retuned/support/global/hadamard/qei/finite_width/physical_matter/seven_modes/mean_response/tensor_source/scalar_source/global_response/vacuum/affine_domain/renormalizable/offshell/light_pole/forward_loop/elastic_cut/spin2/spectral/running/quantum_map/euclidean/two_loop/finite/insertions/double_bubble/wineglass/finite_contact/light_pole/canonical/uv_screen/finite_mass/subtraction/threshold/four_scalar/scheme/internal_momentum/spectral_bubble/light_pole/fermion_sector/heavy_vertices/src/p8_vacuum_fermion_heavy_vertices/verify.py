"""Read-only full nonlocal-vertex fermion-insertion family certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_two_loop_fermion_ledger import verify as parent

from . import audit, bounds, calibration, ownership, parametric, renormalization

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-fermion-heavy-vertices.json"
PARENT_SHA = "0f6a4e972f68a58287a6c3f31faf4b281bf57bc3ed7f0697f083a73141080f5e"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_fermion_heavy_vertices/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen complete primitive fermion-sector ledger changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_137_fully_rebuilt": PARENT_SHA,
        "same_GY14_nongravitational_reference_model": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
        "full_family_bound_not_complete_matching_or_original_P8": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A full heavy-vertex fermion-insertion family gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.138.FULL_HEAVY_VERTEX_FERMION_INSERTION_FAMILY",
        "date": "2026-09-10",
        "status": "COMPLETE_PAIRED_W2_W2_F0_FAMILY_BOUND_IN_DECLARED_LOCAL_PARENT_SUBTRACTION; NOT_COMPLETE_TWO_LOOP_FINITE_MS_MATCHING_HIGHER_LOOP_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/ownership.md",
            "notes/parametric.md",
            "notes/bounds.md",
            "notes/renormalization.md",
            "notes/continuation.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": len(rows),
        "proof_checks": gates,
        "complete_vertex_ownership": serialize(payload(ownership.data())),
        "mass_replacement_and_spectral_bounds": serialize(
            {
                "parameters": payload(parametric.data()),
                "spectral_bounds": payload(bounds.data()),
            }
        ),
        "local_parent_subtraction": serialize(payload(renormalization.data())),
        "actual_enclosure_and_remaining_frontier": serialize(
            {
                "actual": payload(calibration.data()),
                "primitive_frontier": audit.frontier(),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete paired scalar_Phi4_W2_F0 contribution is bounded with both internal insertion positions and the entire nonlocal heavy vertex at each end. One bubble, four triangle and four box products per channel admit a uniform spectral-mass parameter-domain bound. The bubble subtraction is implemented by one fixed set of local full-parent L, G and M counterterms before regulator removal; no independent higher-derivative or tunable b2 contact is used. A rational dyadic logarithm bound gives a uniform amplitude and Cauchy b2 upper bound below 1e-622, relative to the tree coefficient below 1e-22. This replaces the old local-only subset bound for this family, but its finite common-MS/canonical conversion and seven other primitive rows remain unevaluated. It is not the complete two-loop error, higher-loop bound, V contour, finite-gravity G, common-parent B or original P8 closure.",
        "not_established": [
            "Seven remaining primitive fermion-sector rows",
            "Other proper-subgraph counterterms and second-order local references",
            "Finite conversion of this sector to the common MS interaction scheme",
            "Complete second-order canonical fields and fundamental-cubic products",
            "Complete enlarged-model two-loop amplitude or error",
            "Higher-loop errors, V contours, spectrum, cutoff, G, B or original P8 closure",
        ],
        "verification_boundary": "Exact full-Hessian covariance variation and external derivatives, parameter-gap/simplex identities, local parent counterterm algebra, spectral-tail integrals, rational logarithm enclosure and a mutation-checked frontier. Continuum holomorphy and domination are written analytic proofs, not proof-assistant formalization or independent peer review. Native, ordinary, CLI and direct science use unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The full heavy-vertex fermion-insertion family differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.138.FULL_HEAVY_VERTEX_FERMION_INSERTION_FAMILY replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
