"""Read-only complete actual scalar ADM/common-clock reference response."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_heavy_clock_quadratic import verify as clock_input
from p8_vacuum_affine_heavy_curved_state import verify as state_input
from p8_vacuum_affine_heavy_spatial_response import verify as parent

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-heavy-adm-clock-response.json"
)
PARENT_SHA = "b40797cab9d55ff708ec135dc4d7e2c5344bb70480f429ce4c9bba40b986d040"
STATE_SHA = "d4a6081d0482b2168255299cdd0080c90a1716eeffbc70f2645e6e2452e78eda"
CLOCK_SHA = "f7d3b2ddfaf2d8fb2e09dcaf28ca942ddcda8c0ac36ad9e95711b2db9a1e5a04"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_heavy_adm_clock_response/*.py"))
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
        (clock_input, CLOCK_SHA),
    ):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A frozen actual scalar spatial, state or clock input changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_245_entire_actual_scalar_spatial_response_and_original_projection_fully_rebuilt": PARENT_SHA,
        "S6_240_actual_heavy_state_full_covariant_prescription_and_fixed_profile_fully_rebuilt": STATE_SHA,
        "S6_241_actual_nonlinear_common_clock_and_complete_profile_quadratic_fully_rebuilt": CLOCK_SHA,
        "all_frozen_original_and_separate_model_scientific_bytes_unchanged": True,
        "no_old_vector_norm_withdrawn_phase_or_extra_finite_heat_action_transferred": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A complete actual scalar ADM/common-clock gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.246.COMPLETE_ACTUAL_HEAVY_SLE_ADM_AND_COMMON_CLOCK_RESPONSE_WITH_ORDERED_SCALAR_WARD_FULL_METRIC_PROFILE_CONTACTS_AND_ORIGINAL_PROJECTED_MEAN_TAIL_IN_EXPLICIT_WARD_COMPLETED_APPROXIMANT",
        "date": "2026-09-13",
        "status": "COMPLETE_ACTUAL_GAUSSIAN_ADM_COMMON_CLOCK_REFERENCE_RESPONSE_WITH_WARD_COMPLETED_ORIGINAL_PROJECTION; NOT_BARE_CUTOFF_LOCAL_SHAPE_IDENTIFICATION_QUANTUM_INVERSE_NONLINEAR_UV_REGGE_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/geometry.md",
            "notes/ward.md",
            "notes/profile.md",
            "notes/graphs.md",
            "notes/cutoff.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_actual_scalar_metric_and_clock": serialize(
            {
                "entire_metric_scalar_Hamiltonian_and_ordered_Ward_bridge": payload(
                    audit.geometry.data()
                ),
                "full_fixed_profile_and_nonlinear_clock_contacts": payload(
                    audit.profile.data()
                ),
            }
        ),
        "complete_reference_graph_and_restored_projection": serialize(
            payload(audit.estimates.data())
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged actual heavy scalar SLE and full covariant prescription now have their complete first ADM/common-clock response on the reference FLRW slab. The entire spatial response supplies the synchronous Gaussian block only after its full profile QQ term is removed. Both ordered density Ward terms, all nonlinear metric contacts and the complete actual mean remain, and the full fixed ADM profile is restored once. The full nonlinear common-clock contact cancels only after the actual Gaussian and profile means are combined. All ten metric directions and arbitrary external transfers are retained on detector V02 times prepared source U138. The normalized physical ADM and common-clock bounds are below10^-327 and10^-300. The explicitly defined Ward-completed restored original-frequency approximant uses the complete two-leg spatial projection and one-leg mean projection in both ordered Ward terms, retains the nonlinear projected-mean clock tail, and converges on the same graph below10^-200/K. It is not identified with a bare full projected Gaussian response or a proved local/shape counteraction. No compatible quantum inverse, finite inhomogeneous nonlinear feedback, interacting-loop control, quantum gravitational limit, physical UV/Regge or original V/G/B/P8 closure follows.",
        "not_established": [
            "Equality of the Ward-completed computational approximant with a bare full projected Gaussian response or a local/shape counteraction",
            "A compatible quantum constraint inverse, stability or finite inhomogeneous same-state nonlinear feedback",
            "Interacting light/mixed state and omitted-loop control or a full quantum gravitational decoupling limit",
            "A physical momentum cutoff, exact UV S matrix, finite-gravity Regge estimate or all-parent exclusion",
            "A changed state, finite prescription, off-clock parent or cancellation of a projected mean before its clock contact",
            "Original V/G/B or P8 closure",
        ],
        "verification_boundary": "Exact scalar Hamiltonian, full metric/contact and ordered density identities support the written actual-response assembly. Independent literal nonlinear matrix-exponential derivatives, full scalar mode and Hermitian-Gram gauge tangents, nonzero ordered boundary/contact controls and rational whole-graph estimates test all retained terms. Primary scalar retarded metric/contact and conservation context is Hollands-Wald gr-qc/0404074v2 section4.3 equation112 and Theorem5.1; it is not invoked as proof of these numerical bounds or interacting-loop completion. Native/direct/ordinary/CLI retain original SymPy; only full regression uses the separately audited exact-GCD adapter. No frozen scientific bytes are edited. The arguments are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete actual scalar ADM/common-clock report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.246 complete actual heavy ADM/common-clock replay passed; quantum inverse and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
