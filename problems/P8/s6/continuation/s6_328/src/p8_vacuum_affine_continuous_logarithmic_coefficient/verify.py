"""Read-only complete logarithmic coefficient and uniform energy extension."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_known_hard_soft_subtraction import verify as previous
from p8_vacuum_affine_radiative_logarithmic_soft_coefficient import (
    verify as coefficient_input,
)
from p8_vacuum_affine_radiative_state_soft_index import verify as state_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT
    / "certificates/polynomial-vacuum-affine-continuous-logarithmic-coefficient.json"
)
PARENT_SHA = "9d7ab504bf29f8ee50af4d45f6ef7895e88461d09b192ea7f555294971e7a1f8"
COEFFICIENT_SHA = "cb310adf6525b22e6f89eadd2d9badbbfdf1522016831c27d2fedcea9bf806fe"
STATE_SHA = "83fc98b6541ae6fb75c18c342d20c3e4bc03842fe69a3191b3dc8d1a09509df7"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(
            ROOT.glob("src/p8_vacuum_affine_continuous_logarithmic_coefficient/*.py")
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
        (coefficient_input, COEFFICIENT_SHA),
        (state_input, STATE_SHA),
    ):
        if sha(module.REPORT) != digest:
            raise ValueError("A frozen continuous coefficient input changed")
        module.validate_report(
            json.loads(module.REPORT.read_text()), module.build_report()
        )
    return {
        "S6_327_known_hard_soft_and_sharp_log_interference_rebuilt": PARENT_SHA,
        "S6_326_complete_named_radiative_logarithmic_coefficient_rebuilt": COEFFICIENT_SHA,
        "S6_300_original_physical_radiative_recoil_rebuilt": STATE_SHA,
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
        raise ValueError("A continuous coefficient interference proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.328.CONTINUOUS_LOGARITHMIC_COEFFICIENT_AND_UNIFORM_ANGULAR_ENERGY_EXTENSION",
        "date": "2026-09-18",
        "status": "SCOPED_CONTINUOUS_LOG_COEFFICIENT_AND_UNIFORM_ENERGY_EXTENSION; NOT_QUANTUM_ALL_N_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/phase.md",
            "notes/kernels.md",
            "notes/measure.md",
            "notes/calibration.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_continuous_kernels": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_aggregate_phase_and_uniform_measure_extension": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete named S326 coefficient has a unique continuous extension at all additional-soft collinear directions. Exact conservation reorganizes it into jointly continuous single and double angular kernels, including every atom self diagonal, while the aggregate phase depends only on massive recoil. Finite positive angular-energy measures therefore admit a unique extension continuous from weak measure convergence to uniform angular TT norm. Original bounds sharpen to ||C||sup<=293<300 and ||C-C_Born||sup<=23667R<=24000R. This is a kinematic coefficient, not a finite radiative hard-loop amplitude, interacting state, all-N hard sum or P8 closure.",
        "not_established": [
            "Finite radiative hard-loop remainder, unknown hard matching and common regulator",
            "Positive quantum detector probability or interacting quantum state",
            "Summed all-N hard amplitudes and quantum unitarity",
            "Absolute complex Regge, common-parent bounce or original V/G/B/P8 closure",
            "A unique collinear value for an individual null current",
        ],
        "verification_boundary": "Exact generic pair and TT-projector identities, original real/complex finite-state coefficient equality, atom-splitting and phase checks; independent collinear-path and atomic-approximation tests. Written joint-continuity, product-measure and compact finite-net arguments prove the angular-energy extension, not finite samples. Not kernel-formalized. Original SymPy outside FULL.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The continuous coefficient report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.328 continuous logarithmic coefficient and energy extension replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
