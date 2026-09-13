"""Read-only full actual heavy-scalar spatial response with restored original projection."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_heavy_curved_state import verify as state_input
from p8_vacuum_affine_heavy_homogeneous_response import verify as parent
from p8_vacuum_affine_heavy_spatial_uv import verify as uv_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-heavy-spatial-response.json"
PARENT_SHA = "a03aa196b5fb7529ac2c37eee7e79e1e5cec8299ccabbd45e281e0f2f0ce3a32"
STATE_SHA = "d4a6081d0482b2168255299cdd0080c90a1716eeffbc70f2645e6e2452e78eda"
UV_SHA = "1005d258074ca5d2cdd291a5c17010c05710a286f23f5ad1a3064169ec9451e6"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_heavy_spatial_response/*.py"))
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
        (uv_input, UV_SHA),
    ):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A frozen actual scalar state, homogeneous response or UV input changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_244_entire_actual_homogeneous_scalar_response_and_profile_fully_rebuilt": PARENT_SHA,
        "S6_240_actual_heavy_state_and_complete_scalar_prescription_fully_rebuilt": STATE_SHA,
        "S6_243_complete_ordered_scalar_UV_and_same_scheme_finite_difference_fully_rebuilt": UV_SHA,
        "all_frozen_original_and_separate_model_scientific_bytes_unchanged": True,
        "old_vector_multiplicities_finite_weights_and_withdrawn_phases_not_transferred": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A complete actual scalar spatial response gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.245.COMPLETE_ACTUAL_HEAVY_SLE_ALL_TRANSFER_SPATIAL_RESPONSE_WITH_FULL_STATE_TIME_AND_UV_SUBTRACTED_REMAINDERS_DIMENSIONAL_MATCHING_UNCHANGED_PROFILE_AND_RESTORED_ORIGINAL_FREQUENCY_PROJECTION",
        "date": "2026-09-13",
        "status": "COMPLETE_ACTUAL_LINEAR_SPATIAL_RESPONSE_ALL_TRANSFERS_WITH_RESTORED_ORIGINAL_PROJECTION; NOT_FULL_FINITE_CUTOFF_WARD_ADM_CLOCK_ASSEMBLY_QUANTUM_INVERSE_NONLINEAR_UV_REGGE_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/domain.md",
            "notes/state.md",
            "notes/time.md",
            "notes/endpoints.md",
            "notes/dimension.md",
            "notes/assembly.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_scalar_state_joint_domain_and_time_remainders": serialize(
            {
                "whole_joint_domain_and_dimensional_geometry": payload(
                    audit.domain.data()
                ),
                "entire_actual_state_difference_and_six_step_time_remainder": payload(
                    audit.time_remainder.data()
                ),
            }
        ),
        "complete_actual_spatial_remainder_assembly_and_projection": serialize(
            payload(audit.endpoints.data())
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "For the unchanged actual QG2 heavy SLE and complete scalar prescription, the full homogeneous response and profile anchor every spatial transfer. Complete actual-state-minus-W6, sixth-time and UV-subtracted first-five endpoint differences plus the unchanged scalar finite UV difference give the entire linear spatial current. The computational W6 family is not promoted to a physical state. All source jets, trace and tracefree directions, both ordered branches, every endpoint and all initial comparison errors remain. The exact joint complex-radius/time domain controls every full scalar iterate and normalized feature. Complete near/far/low estimates retain all power/log terms; the logarithmic inverse mass is not discarded. Separately normalized scalar geometries and both same-dimension phases justify the full dimensional dominated limit before fixed-invariant finite matching. The complete normalized all-transfer response is below10^-350 on detector L2 times prepared source Z136. The original-frequency two-leg projection applied to the complete unaveraged subtracted integrand has same-graph tail below10^-250/K for K>=2sqrt(n), with all finite local terms and the fixed profile restored unchanged. Angular-zero UV rows are retained before masking. This is not a bare finite-cutoff Ward identity, full ADM/common-clock assembly, quantum inverse, nonlinear inhomogeneous feedback, interacting-loop control, quantum gravitational decoupling, physical UV/Regge or original V/G/B/P8 closure.",
        "not_established": [
            "Full finite-cutoff Ward/local-shape bookkeeping or the complete heavy ADM/common-clock response",
            "A compatible quantum constraints/inverse, stability or nonlinear same-state inhomogeneous bounce",
            "A physical replacement WKB state, changed finite prescription, one-leg cutoff or omitted unaveraged UV terms",
            "Interacting light/mixed state and loop control or quantum gravitational decoupling",
            "A physical momentum cutoff, exact UV S matrix, finite-gravity Regge estimate or all-parent exclusion",
            "Original V/G/B or P8 closure",
        ],
        "verification_boundary": "Exact algebra, full scalar frequency recursions and positive mass-uniform majorants support the written all-transfer arguments. Independent literal sixth-order product rules, high-precision complete retarded integration,39 full unexpanded complex-domain fixtures, exact near/far geometric reconstruction and252 removed-union quadrature cases check signs, coefficients and all momentum regions. Numerical fixtures are not the proof of uniform remainder bounds. Native/direct/ordinary/CLI use original SymPy; only full regression uses the audited exact-GCD adapter. No frozen scientific bytes or prescription are edited. The written arguments are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete actual spatial response report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.245 entire actual heavy spatial response replay passed; full ADM/clock assembly and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
