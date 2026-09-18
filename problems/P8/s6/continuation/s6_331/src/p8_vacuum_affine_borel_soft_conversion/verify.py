"""Read-only Borel soft conversion and conditioned cloud transfer."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_leading_cloud_logarithmic_coefficient import verify as input3
from p8_vacuum_affine_quantitative_soft_cutoff import verify as input0
from p8_vacuum_affine_radiative_angular_finite import verify as input2
from p8_vacuum_affine_radiative_soft_state_transfer import verify as input1

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-borel-soft-conversion.json"
PARENT_SHA = "59468b78c5aa6fa90d162bbad5290d43a2739e4b5310c6d730b64b40bceaca6b"
STATE_SHA = "32396aad07a5c1521a282fe6161c0a8e99ff6e312f55df4937de8f5501e30655"
ANGULAR_SHA = "5e18135d1a49c2846047f2a71c3de00f7c82b84a131d8a6f924890308f4ccd99"
CLOUD_SHA = "34069fadf08a3e44b36cf2ffd7bc909b6dd3fd019fa7ba17c3bfe7d53d567aae"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_borel_soft_conversion/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    for module, digest in (
        (input0, PARENT_SHA),
        (input1, STATE_SHA),
        (input2, ANGULAR_SHA),
        (input3, CLOUD_SHA),
    ):
        if sha(module.REPORT) != digest:
            raise ValueError("A frozen Borel soft conversion input changed")
        module.validate_report(
            json.loads(module.REPORT.read_text()), module.build_report()
        )
    return {
        "S6_330_original_positive_recoil_and_relative_conditioning_rebuilt": PARENT_SHA,
        "S6_325_complete_state_continuity_and_same_state_soft_pairing_rebuilt": STATE_SHA,
        "S6_301_full_fixed_ball_trace_radial_and_phase_prescription_rebuilt": ANGULAR_SHA,
        "S6_329_defined_conditioned_leading_cloud_and_moments_rebuilt": CLOUD_SHA,
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
        raise ValueError("A Borel soft conversion proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.331.POSITIVE_ENERGY_SOFT_CONVERSION_AND_CONDITIONED_CLOUD_TRANSFER",
        "date": "2026-09-18",
        "status": "SCOPED_BOREL_SOFT_CONVERSION_AND_CONDITIONED_LEADING_CLOUD; NOT_HARD_INTERACTING_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/current.md",
            "notes/extension.md",
            "notes/moments.md",
            "notes/regulator.md",
            "notes/cutoff.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_Borel_soft_conversion": serialize(
            {name: payload(packets[name]) for name in names[:3]}
        ),
        "whole_conditioned_cloud_transfer_and_cutoff": serialize(
            {name: payload(packets[name]) for name in names[3:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged additional-soft index, finite conversion and regulated coefficient extend to positive Borel angular-energy measures with explicit positive-increment moduli. On the fixed-Born D4 leading cloud conditioned by total energy, the same-state remaining-energy soft insertion has a total-variation regulator limit<19000x(1+ln1/x)/kappa^2. Its a and Delta marks have physical cutoff-conditioned mean error bounds4000eta/kappa^2 and34000eta(1+ln1/eta)/kappa^2. This does not continue outer hard amplitudes, supply evanescent matching, establish interacting dynamics or close original P8.",
        "not_established": [
            "Positive noninteger-dimensional radiation intensity",
            "Cutoff rate after changing the remaining-energy logarithm inZ",
            "Outer hard dimensional/evanescent matching or finite radiative hard remainder",
            "Full detector probability, state-dependent dynamics or all-N hard sum",
            "Quantum unitarity, absolute Regge, common-parent bounce, UV or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact projected-current algebra and radial majorants, original rational positive-tail calibrations and independent conditional logarithmic integrals support explicit written L2/Borel extension, dominated regulator-limit and relative-conditioning proofs. Numerical samples are not uniform proofs. Only the additional soft factor is dimensionally continued. Not kernel-formalized; original SymPy outside FULL.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The Borel soft conversion report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.331 Borel soft conversion and conditioned cloud transfer replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
