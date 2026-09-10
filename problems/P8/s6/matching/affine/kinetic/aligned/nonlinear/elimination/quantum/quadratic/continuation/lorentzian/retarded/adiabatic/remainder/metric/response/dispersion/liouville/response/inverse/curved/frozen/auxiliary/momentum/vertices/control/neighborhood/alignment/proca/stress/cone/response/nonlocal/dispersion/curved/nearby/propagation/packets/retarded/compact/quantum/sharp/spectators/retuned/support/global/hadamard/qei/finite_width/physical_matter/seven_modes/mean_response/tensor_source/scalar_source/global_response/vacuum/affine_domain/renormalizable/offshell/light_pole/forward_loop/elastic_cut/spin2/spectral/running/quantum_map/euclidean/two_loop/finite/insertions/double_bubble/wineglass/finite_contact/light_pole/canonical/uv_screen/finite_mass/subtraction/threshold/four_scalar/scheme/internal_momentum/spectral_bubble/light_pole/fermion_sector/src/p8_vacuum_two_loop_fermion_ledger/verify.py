"""Read-only complete primitive two-loop fermion-sector ownership ledger."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_fermion_spectral_pole import verify as parent

from . import audit, boundary, counterterms, families, gauge, legendre

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-two-loop-fermion-ledger.json"
PARENT_SHA = "315d9f330c2c92a432cab678646d8b22e847493e2366e1fe343f903a5415902b"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_two_loop_fermion_ledger/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen first spectral quadratic-sector report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_136_fully_rebuilt": PARENT_SHA,
        "same_GY14_nongravitational_reference_model": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
        "complete_ownership_not_complete_estimates_or_original_P8": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A two-loop fermion-sector ownership gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.137.COMPLETE_TWO_LOOP_FERMION_SECTOR_OWNERSHIP",
        "date": "2026-09-10",
        "status": "COMPLETE_PRIMITIVE_FERMION_SECTOR_OWNERSHIP_THROUGH_PHI4; NOT_EVALUATED_TWO_LOOP_ERROR_FINITE_MATCHING_HIGHER_LOOP_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/legendre.md",
            "notes/families.md",
            "notes/gauge.md",
            "notes/counterterms.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": len(rows),
        "proof_checks": gates,
        "loop_and_background_expansion": serialize(
            {"legendre": payload(legendre.data()), "families": payload(families.data())}
        ),
        "gauge_color_flavor_ownership": serialize(payload(gauge.data())),
        "counterterm_and_finite_conversion_ledger": serialize(
            payload(counterterms.data())
        ),
        "exact_remaining_frontier": serialize(payload(boundary.data())),
        "controls": serialize(audit.controls()),
        "verdict": "The complete nongravitational primitive two-loop fermion-sector ownership through four external Phi derivatives is derived by integrating the quadratic fermions, performing the bosonic loop expansion and cancelling reducible source terms with the Legendre transform. The noncommuting scalar expansion has six rows and the gauge contraction has three, with all fourteen vacuum flavors and distinct active scalar/gauge color factors retained. Separate bosonic/fermionic one-loop counterterm insertions, second-order references, regulator-sensitive finite products and canonical/fundamental-cubic conversion terms are explicit tasks; assigned proper counterterms may not be counted twice. The quantitative frontier remains incomplete: one paired quadratic row is bounded, one quartic row has only a local-vertex subset bounded, and seven primitive rows plus the other counterterm/conversion tasks are unevaluated. A complete catalog is not a complete two-loop amplitude or pole/error budget, later-loop bound, V contour, finite-gravity G, common-parent B or original P8 closure.",
        "not_established": [
            "Seven primitive rows and the remaining nonlocal heavy-vertex quartic subset",
            "All remaining proper-subgraph and second-order reference coefficients",
            "Complete finite MS interaction and canonical-field conversion",
            "Complete enlarged-model two-loop amplitude, pole or error budget",
            "Higher-loop/truncation errors, V contours, gauge spectrum and cutoff",
            "Finite-gravity Regge/IR Delta for G and controlled common-parent B",
            "Original P8 closure or exclusion of a complete old ladder row",
        ],
        "verification_boundary": "Exact Gaussian moments and Legendre algebra, noncommuting matrix/background expansion, literal SU(3) color/flavor identities, proper-counterterm partition rules, finite regulator products and a mutation-checked quantitative frontier. The arbitrary-index/regulated functional argument is a written proof cross-checked against primary literature, not proof-assistant formalization, evaluated continuum errors or an all-orders QFT theorem. Native and ordinary replay uses unmodified SymPy with interpreter-only runtime allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The two-loop fermion-sector ledger differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.137.COMPLETE_TWO_LOOP_FERMION_SECTOR_OWNERSHIP replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
