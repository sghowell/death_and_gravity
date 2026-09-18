"""Read-only positive-energy coefficient and quantitative infrared cutoff."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_continuous_logarithmic_coefficient import (
    verify as coefficient_input,
)
from p8_vacuum_affine_leading_cloud_logarithmic_coefficient import verify as previous
from p8_vacuum_affine_radiative_logarithmic_soft_coefficient import verify as pair_input
from p8_vacuum_affine_radiative_state_soft_index import verify as recoil_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-quantitative-soft-cutoff.json"
PARENT_SHA = "34069fadf08a3e44b36cf2ffd7bc909b6dd3fd019fa7ba17c3bfe7d53d567aae"
COEFFICIENT_SHA = "bac38028285504de71100c03dd1ed12667b0d0f89e6f79130405be6f05d33fb8"
PAIRS_SHA = "cb310adf6525b22e6f89eadd2d9badbbfdf1522016831c27d2fedcea9bf806fe"
RECOIL_SHA = "83fc98b6541ae6fb75c18c342d20c3e4bc03842fe69a3191b3dc8d1a09509df7"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_quantitative_soft_cutoff/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    for module, digest in (
        (previous, PARENT_SHA),
        (coefficient_input, COEFFICIENT_SHA),
        (pair_input, PAIRS_SHA),
        (recoil_input, RECOIL_SHA),
    ):
        if sha(module.REPORT) != digest:
            raise ValueError("A frozen quantitative cutoff input changed")
        module.validate_report(
            json.loads(module.REPORT.read_text()), module.build_report()
        )
    return {
        "S6_329_conditioned_leading_cloud_and_moments_rebuilt": PARENT_SHA,
        "S6_328_complete_continuous_coefficient_and_aggregate_phase_rebuilt": COEFFICIENT_SHA,
        "S6_326_uniform_original_pair_and_phase_bounds_rebuilt": PAIRS_SHA,
        "S6_300_exact_original_positive_radiative_recoil_rebuilt": RECOIL_SHA,
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
        raise ValueError("A quantitative cutoff interference proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.330.POSITIVE_ENERGY_LIPSCHITZ_COEFFICIENT_AND_QUANTITATIVE_INFRARED_CUTOFF",
        "date": "2026-09-18",
        "status": "SCOPED_POSITIVE_ENERGY_COEFFICIENT_AND_QUANTITATIVE_LEADING_CUTOFF; NOT_INTERACTING_ALL_N_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/increment.md",
            "notes/cutoff.md",
            "notes/conditioning.md",
            "notes/calibration.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_positive_energy_increment": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_quantitative_conditioned_infrared_cutoff": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete coefficient is uniformly Lipschitz under positive added angular energy at fixed original E,u: ||C_sigma-C_tau||sup<=40000*mass(sigma-tau). In the same fixed-Born D4 leading cloud, exact distinguished-emission tail moments and a relative total-cut normalization estimate yield physical mean cutoff error<70000*eta/kappa^(5/2) and second-moment cutoff error<2000000000*x*eta/kappa^4. These bounds remain uniform for0<eta<=x<=1/8 without expanding the tiny calorimetric normalization or erasing the complex phase. They do not supply signed-measure perturbations, marked dimensional matching, hard radiative amplitudes, interacting states or original P8 closure.",
        "not_established": [
            "Arbitrary signed perturbations or changing hard E,u",
            "Marked dimensional/evanescent conversion or common-regulator detector matching",
            "Finite radiative hard remainder, interacting state or all-N hard summation",
            "Positive full detector rate, quantum unitarity or absolute complex Regge",
            "Common-parent bounce, UV completion or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact original recoil rewrite, coefficient-difference algebra, derivative budgets and Poisson/Bose tail identities, with independent original positive-tail calibrations and conditional-moment tests. Written compact-domain derivative and relative-conditioning proofs establish uniformity; finite samples do not. Only the positive-energy increment and fixed-Born D4 marked reference are covered. Not kernel-formalized. Original SymPy outside FULL.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The quantitative cutoff report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.330 positive-energy and quantitative infrared cutoff replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
