"""Read-only remaining-energy common-space soft-cutoff report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_borel_soft_conversion import verify as input0
from p8_vacuum_affine_leading_cloud_logarithmic_coefficient import verify as input3
from p8_vacuum_affine_quantitative_soft_cutoff import verify as input1
from p8_vacuum_affine_radiative_soft_state_transfer import verify as input2

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-remaining-log-cutoff.json"
PARENT_SHA = "642539fc5235de3a4a79099a5f105785cd89535e14ec3a46636f5957cfe3bea8"
CUTOFF_SHA = "59468b78c5aa6fa90d162bbad5290d43a2739e4b5310c6d730b64b40bceaca6b"
STATE_SHA = "32396aad07a5c1521a282fe6161c0a8e99ff6e312f55df4937de8f5501e30655"
CLOUD_SHA = "34069fadf08a3e44b36cf2ffd7bc909b6dd3fd019fa7ba17c3bfe7d53d567aae"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_remaining_log_cutoff/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    for module, digest in (
        (input0, PARENT_SHA),
        (input1, CUTOFF_SHA),
        (input2, STATE_SHA),
        (input3, CLOUD_SHA),
    ):
        if sha(module.REPORT) != digest:
            raise ValueError("A frozen remaining-log cutoff input changed")
        module.validate_report(
            json.loads(module.REPORT.read_text()), module.build_report()
        )
    return {
        "S6_331_positive_Borel_moduli_and_same_state_insertion_rebuilt": PARENT_SHA,
        "S6_330_relative_cut_conditioning_and_renewal_density_rebuilt": CUTOFF_SHA,
        "S6_325_actual_Born_moduli_and_additional_soft_prescription_rebuilt": STATE_SHA,
        "S6_329_defined_leading_reference_and_exact_energy_moments_rebuilt": CLOUD_SHA,
        "S279_remains_rejected_and_archived": True,
        "S275_S276_S277_and_scoped_P8a_unchanged": True,
        "all6_historical_physical_qualifications_unchanged": True,
        "all_original_primitive_and_matching_frontiers_retained": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A remaining-log cutoff proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.332.REMAINING_ENERGY_LOG_COMMON_SPACE_SOFT_CUTOFF",
        "date": "2026-09-18",
        "status": "SCOPED_QUANTITATIVE_COMMON_SPACE_SOFT_MARK_CUTOFF; NOT_HARD_INTERACTING_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/density.md",
            "notes/tail.md",
            "notes/endpoint.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_density": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_same_event_and_excluded_log": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The same limiting additional-soft insertion, including its actual remaining-energy logarithm and exact changed conditioning, has common-space signed-density L1 cutoff error <=63000*a*eta*(1+ln(1/eta))^2/kappa, physically <51000*eta*(1+ln(1/eta))^2/kappa^2. This completes the specific limiting-log cutoff estimate left open in S331. It does not assert TV convergence of emission-configuration laws, joint regulator control, hard matching, interacting dynamics or original P8 closure.",
        "not_established": [
            "TV convergence of finite and infinite emission-configuration laws",
            "Simultaneous quantitative epsilon and eta regulator estimate",
            "Outer hard dimensional/evanescent matching or finite radiative hard remainder",
            "Full detector probability, state-dependent dynamics or all-N hard sum",
            "Quantum unitarity, absolute Regge, common-parent bounce, UV or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact original-source modulus checks, density normalization and logarithmic calculus support written distinguished-emission, renewal-density and common-space L1 proofs. Independent quadratures and finite-Poisson model controls are calibrations, not uniform proofs or replacements of the original current. Not kernel-formalized; original SymPy outside FULL.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The remaining-log cutoff report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.332 remaining-energy common-space soft cutoff replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
