"""Read-only entire actual heavy homogeneous response and full fixed-profile certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_adiabatic_action import verify as action_input
from p8_vacuum_affine_heavy_curved_state import verify as state_input
from p8_vacuum_affine_heavy_spatial_uv import verify as parent

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-heavy-homogeneous-response.json"
)
PARENT_SHA = "1005d258074ca5d2cdd291a5c17010c05710a286f23f5ad1a3064169ec9451e6"
STATE_SHA = "d4a6081d0482b2168255299cdd0080c90a1716eeffbc70f2645e6e2452e78eda"
ACTION_SHA = "98888c7a7fdcfe4fee142833a4220c62192b2834f736625192c8534201dca278"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_heavy_homogeneous_response/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    for prior, expected in (
        (parent, PARENT_SHA),
        (state_input, STATE_SHA),
        (action_input, ACTION_SHA),
    ):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A frozen scalar state, UV or universal variational input changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_243_complete_scalar_spatial_UV_input_fully_rebuilt": PARENT_SHA,
        "S6_240_actual_state_full_scalar_prescription_and_fixed_profile_fully_rebuilt": STATE_SHA,
        "S6_193_universal_single_oscillator_current_action_identity_only_fully_rebuilt": ACTION_SHA,
        "all_frozen_original_and_separate_model_scientific_bytes_unchanged": True,
        "old_vector_multiplicities_finite_weights_and_withdrawn_phases_not_transferred": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A complete actual homogeneous scalar response gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.244.COMPLETE_ACTUAL_HEAVY_SLE_HOMOGENEOUS_SPATIAL_RESPONSE_WITH_FULL_PREPARED_COVARIANCE_AND_CONTACTS_MASS_UNIFORM_ENTIRE_REMAINDERS_GENERAL_COVARIANT_MATCHING_AND_UNCHANGED_FIXED_PROFILE",
        "date": "2026-09-13",
        "status": "COMPLETE_ACTUAL_HOMOGENEOUS_SPATIAL_RESPONSE_WITH_FULL_SCALAR_PRESCRIPTION_AND_FIXED_PROFILE; NOT_NONZERO_TRANSFER_ADM_CLOCK_ASSEMBLY_QUANTUM_INVERSE_NONLINEAR_UV_REGGE_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/evolution.md",
            "notes/preparation.md",
            "notes/bounds.md",
            "notes/matching.md",
            "notes/local.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_actual_evolution_preparation_and_mass_uniform_remainders": serialize(
            {
                "whole_scalar_evolution_and_current": payload(audit.evolution.data()),
                "entire_uniform_prepared_state_and_subtraction_bounds": payload(
                    audit.estimates.data()
                ),
            }
        ),
        "complete_covariant_matching_finite_profile_and_homogeneous_graph": serialize(
            payload(audit.renormalization.data())
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "For the unchanged actual QG2 heavy SLE, the full homogeneous scalar Hamiltonian and time-dependent canonical frame give the entire Gaussian covariance and current. Both first and second prepared source derivatives retain every metric contact and the quadratic first-error term. The full initial state-selection, exact-comparison and W6-to-finite-reference mismatch is bounded rather than reset. Exact mixed-jet majorants and a mass-independent auxiliary marker circle bound the entire actual-state-minus-adiabatic current and its first two source derivatives at all internal momenta throughout a finite prescribed C12 history tube. General-dimensional Cartesian/radial action matching, with arbitrary noncommuting extrinsic-curvature jets and all weighted boundaries, yields exactly the original scalar finite prescription. Its whole local Euler current and the complete unchanged fixed profile are included once. The complete normalized homogeneous first spatial response is below10^-350 on detector time L2 times prepared source H13. The same original-frequency projection has a same-graph tail below10^-748/K squared for K>=2sqrt(n); all finite local/profile terms remain full and no finite-K matched mean is set to zero. This is the entire zero-transfer anchor, not only the S242 state-selection difference or a local heat expansion. Nonzero-transfer remainders, full ADM/common-clock Ward assembly, quantum constraints/inverse, light/mixed loops, nonlinear same-state bounce, quantum gravitational decoupling, physical UV, Regge and original V/G/B/P8 remain open.",
        "not_established": [
            "The full nonzero-transfer UV-subtracted response or its complete original-regulator remainder",
            "The full heavy ADM/common-clock Ward assembly, quantum constraints/inverse or stable nonlinear self-consistent bounce",
            "A replacement physical WKB state, live reminimization, omitted initial mismatch or changed finite prescription",
            "Smallness of a same-space inverse, interacting light/mixed loops or quantum gravitational decoupling",
            "A physical momentum cutoff, exact UV S matrix, finite-gravity Regge estimate or all-parent exclusion",
            "Original V/G/B or P8 closure",
        ],
        "verification_boundary": "Exact generic algebra and positive rational majorants support the written all-momentum arguments. Independent literal Christoffel/Cartesian moments, high-precision full W6 initial tails, direct scalar covariance evolution and finite-difference current derivatives check the derivation with nonzero contacts. These diagnostics are not the all-momentum proof. Native/direct/ordinary/CLI use original SymPy; only full regression uses the audited exact-GCD adapter. No frozen scientific bytes or prescription are edited. The written proofs are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete actual homogeneous response report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.244 entire actual homogeneous heavy response replay passed; nonzero transfer and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
