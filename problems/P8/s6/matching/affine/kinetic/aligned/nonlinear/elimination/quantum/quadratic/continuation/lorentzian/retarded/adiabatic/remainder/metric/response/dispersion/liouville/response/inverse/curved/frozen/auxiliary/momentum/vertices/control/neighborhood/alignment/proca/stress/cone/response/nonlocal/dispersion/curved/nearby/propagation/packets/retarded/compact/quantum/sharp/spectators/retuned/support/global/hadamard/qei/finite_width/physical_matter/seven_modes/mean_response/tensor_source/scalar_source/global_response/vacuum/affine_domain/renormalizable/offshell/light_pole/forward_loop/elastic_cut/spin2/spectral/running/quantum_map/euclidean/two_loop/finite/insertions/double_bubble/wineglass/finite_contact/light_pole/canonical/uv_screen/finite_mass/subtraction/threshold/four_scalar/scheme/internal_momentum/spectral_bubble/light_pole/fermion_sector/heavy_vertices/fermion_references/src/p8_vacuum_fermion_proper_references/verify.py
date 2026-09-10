"""Read-only proper one-loop fermion local-reference certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_fermion_heavy_vertices import verify as parent

from . import anchors, audit, bounds, calibration, conversion, dirac

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-fermion-proper-references.json"
PARENT_SHA = "478b3e99cd4badc723c854a378d17adbd58dd0e9b6252aec1d358140e10805eb"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_fermion_proper_references/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen full heavy-vertex insertion-family report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_138_fully_rebuilt": PARENT_SHA,
        "same_GY14_nongravitational_reference_model": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
        "proper_references_not_complete_matching_or_original_P8": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A proper fermion local-reference gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.139.FERMION_PROPER_LOCAL_REFERENCE_MATCHING",
        "date": "2026-09-10",
        "status": "ONE_LOOP_FERMION_PROPER_LOCAL_REFERENCES_AND_FINITE_CANONICAL_DICTIONARY; NOT_COMPLETE_TWO_LOOP_PRIMITIVES_OUTER_REFERENCE_CONVERSION_HIGHER_LOOP_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/dirac.md",
            "notes/anchors.md",
            "notes/bounds.md",
            "notes/counterterms.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "fermion_inverse_and_dimensional_anchors": serialize(
            {"dirac": payload(dirac.data()), "anchors": payload(anchors.data())}
        ),
        "proper_counterterms_and_field_conversion": serialize(
            payload(conversion.data())
        ),
        "exact_reference_enclosures": serialize(payload(bounds.data())),
        "actual_common_reference_and_primitive_frontier": serialize(
            {
                "actual": payload(calibration.data()),
                "primitive_frontier": audit.previous.frontier(),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The scalar and Feynman-gauge one-loop proper fermion mass, kinetic and scalar-Yukawa local references are derived with the open-line color factor, complete d-dimensional numerator and finite epsilon-times-pole terms. Fixed-scale background differentiation gives the zero-momentum Yukawa vertex; its UV counterterm, including both fermion fields and the scalar field, reproduces the independent frozen matrix beta function. Explicit finite MS-to-zero-momentum canonical ratios retain the same complete physical Phi normalization. Exact rational bounds give relative selected mass and Yukawa shifts below 1e-204. Active opposite Yukawa signs and inert gauge-only flavors are distinguished. These are local proper-reference data, not newly evaluated two-loop primitive integrals, the previous outer-reference finite conversion, a charged-fermion physical spectrum, a complete matching error, V, G, B or original P8 closure.",
        "not_established": [
            "Seven remaining primitive fermion-sector integrals and their complete forests",
            "Complete order-two counterterm-insertion integrals and finite references",
            "S6.138 paired-zero-bubble finite conversion to the common interaction scheme",
            "Complete order-two canonical matching, pole or amplitude error",
            "Exact charged-fermion spectrum or gauge-independent off-shell coefficients",
            "Higher-loop errors, V contours, cutoff, G, B or original P8 closure",
        ],
        "verification_boundary": "Literal Euclidean Clifford/open-color identities and Gaussian inverse-resolvent signs, exact scalar parameter moments, dimensional finite products, fixed-scale mass differentiation, independently cross-checked Yukawa UV running, finite canonical ratios and rational bounds. The continuum local-subgraph derivation is a written proof, not proof-assistant formalization or evaluated two-loop continuum errors. Native, ordinary, CLI and direct science use unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The proper fermion local-reference differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.139.FERMION_PROPER_LOCAL_REFERENCE_MATCHING replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
