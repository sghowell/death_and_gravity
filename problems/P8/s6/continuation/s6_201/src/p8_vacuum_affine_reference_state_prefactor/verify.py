"""Read-only curved initial-state prefactor removal and full response correction."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_flat_spatial_conversion import verify as parent
from p8_vacuum_affine_proca_gaussian import verify as gaussian

from . import audit, contact, limit, memory, prefactor

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-reference-state-prefactor.json"
)
PARENT_SHA = "739ab8954fc755149312818d2c6db8e5d38d943a65cb085160a337891bf30507"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_reference_state_prefactor/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    for previous, expected in ((parent, PARENT_SHA), (gaussian, GAUSSIAN_SHA)):
        if sha(previous.REPORT) != expected:
            raise ValueError(
                "The source-pinned actual affine Gaussian state chain changed"
            )
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_200_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A curved reference initial-state prefactor gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.201.ACTUAL_CURVED_PROCA_INITIAL_STATE_PREFACTOR_REMOVAL_WITH_FULL_MEMORY_CONTACT_TAIL_AND_UNIT_W8_REFERENCE_DECOMPOSITION",
        "date": "2026-09-11",
        "status": "ACTUAL_CURVED_PROCA_INITIAL_STATE_PREFACTOR_REMOVAL_WITH_FULL_MEMORY_CONTACT_TAIL_AND_UNIT_W8_REFERENCE_DECOMPOSITION; NOT_JOINT_SPATIAL_ANALYTICITY_COVARIANT_MATCHING_FULL_INVERSE_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/state.md",
            "notes/memory.md",
            "notes/contact.md",
            "notes/tail.md",
            "notes/reference.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_initial_occupation_and_retarded_memory_correction": serialize(
            {"prefactor": payload(prefactor.data()), "memory": payload(memory.data())}
        ),
        "full_contact_regulator_limit_and_curved_reference_decomposition": serialize(
            {"contact": payload(contact.data()), "limit": payload(limit.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged actual CD prepared state has constant initial phases that cancel exactly in the complete quadratic reference response. All occupation weights, including the quadratic two-leg product, give a full overlapping-time memory and noncommuting metric-contact correction below2e12 M[D]M[Gamma], with the original distinct-regulator error below1e17 M[D]M[Gamma]/K. No momentum derivative or analytic assumption on the Borel state prefactor is used. Combining the existing actual-state remainder gives a unit-W8 comparison below3e23 and error2e27/K. Repeating the full six-step identity with the same amplitude upper bounds retains a known finite curved piece below2e48 and error2e52/K, with canonical displays8e-752 and8e-748/K. Unit W8 has the canonical Wronskian but is not an exact bisolution or a replacement physical state. The remaining full contact and first five equal-time kernels still need joint spatial estimates and fixed covariant matching. This does not establish a full inverse, nonlinear background, stability, physical cutoff, remaining parent loops or finite-gravity IR/Regge. Original V/G/B and P8 remain open.",
        "not_established": [
            "Joint spatial momentum analyticity or the full reference spatial Taylor/subtraction expansion",
            "Original fixed covariant matching of the complete reference contact and first five equal-time kernels",
            "A unit-W8 exact bisolution or new physical state, full mixed inverse, interacting background or stability",
            "Full parent/cutoff/scattering/IR/Regge matching or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact complete initial-phase cancellation, both-leg occupation weights, ten-field contact, massive radial constants, current telescoping and canonical factors. Independent complex readouts, noncommuting contact matrices, Wronskian and non-bisolution fixtures, radial integrals and full removed-union tails supplement written continuum estimates. Native, direct science, ordinary and CLI use original SymPy; only full regression uses the audited exact GCD adapter. The continuous bounds and reference-decomposition argument are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The curved reference initial-state prefactor report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.201.ACTUAL_CURVED_PROCA_INITIAL_STATE_PREFACTOR_REMOVAL_WITH_FULL_MEMORY_CONTACT_TAIL_AND_UNIT_W8_REFERENCE_DECOMPOSITION replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
