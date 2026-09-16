"""Read-only full leading-soft calorimetric resummation certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_physical_virtual_soft_pairing import verify as soft_input
from p8_vacuum_affine_unequal_mass_regge_residue import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-calorimetric-soft-resummation.json"
)
PARENT_SHA = "322c3cfc4ed1bcb9c209d54ead72327b694b2e259a5525d2f0d067ad675c83aa"
SOFT_SHA = "40c04cadfa1a1c542c16eed48f6df12bcfe0605ad69b39dbb70fcefaefe6b3ed"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_calorimetric_soft_resummation/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen unequal-mass Regge parent report changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(soft_input.REPORT) != SOFT_SHA:
        raise ValueError("The frozen complete physical soft conversion changed")
    soft_input.validate_report(
        json.loads(soft_input.REPORT.read_text()), soft_input.build_report()
    )
    return {
        "S6_298_unequal_mass_Regge_ancestry_rebuilt": PARENT_SHA,
        "S6_296_complete_D_angular_conversion_rebuilt": SOFT_SHA,
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
        raise ValueError("A complete leading-soft calorimetric proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.299.COMPLETE_LEADING_SOFT_CALORIMETRIC_SUM_FINITE_ANGULAR_CONVERSION_AND_ORDERED_SCALING",
        "date": "2026-09-15",
        "status": "SCOPED_COMPLETE_ALL_LEADING_SOFT_CALORIMETRIC_SUM_AND_ORDERED_LIMITS; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/poisson.md",
            "notes/gamma.md",
            "notes/limits.md",
            "notes/matching.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_total_energy_sum": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_Gamma_conversion_and_detector_limits": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The total-energy simplex sums every leading-soft graviton multiplicity. A positive auxiliary Poisson representation rigorously yields exp(Delta)*exp(-EulerGamma*a)/Gamma(1+a)*(E/nu)^a after removing the dimensional regulator. The finite angular conversion is exactly the inherited S296 value. The whole Gamma-product inequality bounds the original conversion remainder beyond1+Delta below2*10^-1596, uniformly in the unexpanded resolution power. Ordered and joint-regulator detector scalings are explicitly distinguished. Uniform full-amplitude radiation and hard/Regge/matching errors remain unproved; original V/G/B/P8 remain OPEN.",
        "not_established": [
            "Uniform all-multiplicity full-radiation-minus-leading-soft control",
            "Gravity-Born recoil and all omitted soft or hard loop bounds",
            "Finite gravitational matching and full physical analytic-amplitude relation",
            "High-energy complex-energy/Regge remainder control",
            "Original state/domain/measure/bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "The written Dirichlet/Poisson proof controls the complete defined leading-soft positive sum; it is not an exchange of its approximation with the full infinite-multiplicity radiation problem. Exact algebra checks normalization, all original inputs, Gamma inequalities and limit scalings. Independent finite-regulator sums and Laplace inversion calibrate the whole functions, not a substitute for the written bounds. Native/direct/ordinary/CLI use original SymPy; the adapter remains FULL-only. All frozen ancestors, historical qualifications and original closure statuses remain unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete leading-soft calorimetric report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.299 complete leading-soft calorimetric replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
