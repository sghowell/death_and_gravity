"""Read-only complete minimally coupled scalar spatial UV and finite-difference certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_heavy_curved_state import verify as state_input
from p8_vacuum_affine_heavy_state_response import verify as parent
from p8_vacuum_affine_ordered_scalar_symbol import verify as geometry_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-heavy-spatial-uv.json"
PARENT_SHA = "a9a26c7e913edd1e6fb8a4455994066c884bc4fe7979b168fe7871fcef6059b4"
STATE_SHA = "d4a6081d0482b2168255299cdd0080c90a1716eeffbc70f2645e6e2452e78eda"
GEOMETRY_SHA = "70e8538425cfab3035cea3c4c9d3a5da63797182625f79811bbc9ad2cd9ed815"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_heavy_spatial_uv/*.py"))
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
        (geometry_input, GEOMETRY_SHA),
    ):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A frozen QG2 state, response or generic geometry input changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_242_complete_state_selection_response_and_fixed_profile_fully_rebuilt": PARENT_SHA,
        "S6_240_actual_state_and_entire_covariant_scalar_prescription_fully_rebuilt": STATE_SHA,
        "S6_216_corrected_generic_full_six_invariant_geometric_Hessians_fully_rebuilt": GEOMETRY_SHA,
        "all_frozen_original_and_separate_model_scientific_bytes_unchanged": True,
        "old_vector_physical_weights_and_withdrawn_phase_formulas_not_reendorsed": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "A complete scalar spatial UV or finite-difference gate failed"
        )
    return {
        "schema": 1,
        "claim": "P8-S6.243.COMPLETE_MINIMAL_SCALAR_ORDERED_SPATIAL_UV_SYMBOL_WITH_ALL_SIX_INVARIANTS_FULL_COVARIANT_POLE_AND_SAME_SCHEME_DIMENSIONAL_FINITE_DIFFERENCE_AND_LOCAL_GRAPH_BOUND",
        "date": "2026-09-13",
        "status": "COMPLETE_LOCAL_SPATIAL_SCALAR_UV_INPUT_AND_SAME_SCHEME_FINITE_DIFFERENCE; NOT_COMPLETE_UV_SUBTRACTED_COMPARISON_RESPONSE_HOMOGENEOUS_ANCHOR_QUANTUM_INVERSE_NONLINEAR_UV_REGGE_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/vertex.md",
            "notes/ordering.md",
            "notes/dimension.md",
            "notes/finite.md",
            "notes/bounds.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_ordered_scalar_UV_jet_and_dimensional_input": serialize(
            {
                "entire_compiled_jet": audit.complete_jet_payload(),
                "full_vertex_ordering_and_normalization": payload(
                    audit.structure.data()
                ),
            }
        ),
        "complete_covariant_finite_spatial_difference_and_bound": serialize(
            payload(audit.finite_input.finite())
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "For the unchanged exact QG2 heavy state and its actual minimal-scalar covariant prescription, the full spatial Hamiltonian vertex is derived before UV extraction. Both ordered trace/gradient cross products, all60 endpoint/source slots and140 coefficients are retained with the physical annihilation phase. The complete35 radial rows include all time, mass, trace, gradient and angular contributions. The arbitrary-dimensional transverse average keeps all six physical invariants fixed before differentiation. The entire scalar logarithmic symbol matches the full covariant scalar pole action in every invariant; other apparent spatial power/odd grades vanish only after their complete ordered sum. Exact dimensional sphere and MSbar normalization plus the full volume/Euler evanescent counteraction give the actual same-scheme finite UV difference, not a fitted target or the old vector finite weights. All ordered local Green identities hold. The full eighteen-coefficient rational/logarithmic sum proves a normalized local bound below10^-600 on detector L2 times source Z24 for all transfers on the clock slab. This is the complete LOCAL spatial UV input, not the complete UV-subtracted exact-comparison kernel, homogeneous anchor or full heavy ADM/clock response. The S241 finite heat Hessian is not added again as another counterterm. No quantum inverse, stability, interacting loop control, nonlinear bounce, quantum gravitational decoupling, physical UV, finite-gravity Regge or original V/G/B/P8 closure follows.",
        "not_established": [
            "The full UV-subtracted exact-comparison response, its homogeneous anchor, regulator tails or full-dimensional dominated limit",
            "The complete heavy ADM/common-clock response, quantum constraints/inverse or nonlinear same-state bounce",
            "A replacement physical WKB state, live reminimization or changed finite prescription",
            "Smallness of a same-space inverse, interacting light/mixed loops or quantum gravitational decoupling",
            "A physical cutoff, exact UV S matrix, finite-gravity Regge estimate or all-parent exclusion",
            "Original V/G/B or P8 closure",
        ],
        "verification_boundary": "Exact generic algebra derives the complete UV input and every finite coefficient. Independent Cartesian spherical moments, literal scalar Hamiltonian features, full numerical six-iterate time jets and dimensional Gamma normalization cross-check it; numerical fixtures do not provide all-momentum remainder estimates. Native/direct/ordinary/CLI use original SymPy; only full regression uses the audited exact-GCD adapter. No frozen scientific bytes or prescription are edited. The written proofs are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete scalar spatial UV report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.243 complete scalar spatial UV replay passed; full comparison response and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
