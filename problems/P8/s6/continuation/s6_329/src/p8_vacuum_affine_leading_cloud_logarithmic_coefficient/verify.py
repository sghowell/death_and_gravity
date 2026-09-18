"""Read-only conditioned leading-cloud logarithmic coefficient report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_calorimetric_soft_resummation import verify as leading_input
from p8_vacuum_affine_continuous_logarithmic_coefficient import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT
    / "certificates/polynomial-vacuum-affine-leading-cloud-logarithmic-coefficient.json"
)
PARENT_SHA = "bac38028285504de71100c03dd1ed12667b0d0f89e6f79130405be6f05d33fb8"
LEADING_SHA = "604445d8362c3dd9cec6e456bab0743e0611aaac2e60794e64f99be41cd63879"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(
            ROOT.glob("src/p8_vacuum_affine_leading_cloud_logarithmic_coefficient/*.py")
        )
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
        (leading_input, LEADING_SHA),
    ):
        if sha(module.REPORT) != digest:
            raise ValueError("A frozen leading-cloud coefficient input changed")
        module.validate_report(
            json.loads(module.REPORT.read_text()), module.build_report()
        )
    return {
        "S6_328_continuous_complete_coefficient_and_energy_extension_rebuilt": PARENT_SHA,
        "S6_299_full_calorimetric_leading_reference_and_original_index_bound_rebuilt": LEADING_SHA,
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
        raise ValueError("A leading-cloud coefficient interference proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.329.CONDITIONED_LEADING_CLOUD_LOGARITHMIC_COEFFICIENT_MOMENTS",
        "date": "2026-09-18",
        "status": "SCOPED_CONDITIONED_LEADING_CLOUD_COEFFICIENT; NOT_INTERACTING_ALL_N_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/cloud.md",
            "notes/conditioning.md",
            "notes/moments.md",
            "notes/calibration.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_defined_leading_cloud": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_conditioned_complete_coefficient_moments": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete continuous logarithmic coefficient is a well-defined mark on the fixed-Born D4 leading soft cloud conditioned on total energy at most x. Its infrared-cutoff limit exists in every finite Lp of the uniform angular TT norm, including the infinite-particle finite-energy limit. Exact conditional energy moments give physical coefficient mean change<19000*x/kappa^(5/2) and second moment<230000000*x^2/kappa^4. Full complex phases and the unexpanded calorimetric normalization are kept. This is not an interacting state, a marked dimensional-scheme conversion, all-N hard summation or original P8 closure.",
        "not_established": [
            "Full radiative hard amplitude and common-regulator detector matching",
            "Dimensional or evanescent conversion for the marked insertion",
            "State-dependent branching law or interacting quantum state",
            "Summed all-N hard amplitudes, positive full detector probability or unitarity",
            "Absolute complex Regge, common-parent bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact original leading CDF, Levy exponent, cutoff and conditioned-moment algebra and original-parameter bounds. The written independent-band probability construction, nested-event conditioning and S328 uniform-continuity argument prove the infinite-cloud and marked Lp limits. Finite tests do not construct an interacting quantum state; no Monte Carlo inference is used for the extremely small physical index. Not kernel-formalized. Original SymPy outside FULL.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The leading-cloud coefficient report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.329 conditioned leading-cloud coefficient replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
