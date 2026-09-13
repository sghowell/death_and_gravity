"""Read-only actual current conditional two-mass coupled inverse report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_band_coupled_inverse import verify as generic_input
from p8_vacuum_affine_heavy_clock_quadratic import verify as current_input
from p8_vacuum_affine_two_mass_reference import verify as parent

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-two-mass-coupled-inverse.json"
)
PARENT_SHA = "dce8b1593e15778f9b5ecb2ea224e0cb9b834d2e23a6d1c6fad126c358a8afd0"
CURRENT_SHA = "f7d3b2ddfaf2d8fb2e09dcaf28ca942ddcda8c0ac36ad9e95711b2db9a1e5a04"
GENERIC_SHA = "076648eee2112f5b86302a072e11265350e155a92877f5b68b5007ab0cf09432"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_two_mass_coupled_inverse/*.py"))
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
        (current_input, CURRENT_SHA),
        (generic_input, GENERIC_SHA),
    ):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A frozen actual state, current coefficient or complete reference input changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_247_complete_same_scheme_two_mass_reference_and_full_response_ancestry_fully_rebuilt": PARENT_SHA,
        "S6_241_whole_current_QG2_classical_coefficients_and_profile_fully_rebuilt": CURRENT_SHA,
        "S6_229_generic_coordinate_and_Proca_remainder_input_fully_rebuilt": GENERIC_SHA,
        "all_frozen_original_and_separate_model_scientific_bytes_unchanged": True,
        "old_QG1_inverse_status_not_transferred_to_the_new_actual_operator": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("An actual two-mass coupled inverse gate failed")
    packets = audit.packets()
    return {
        "schema": 1,
        "claim": "P8-S6.248.ACTUAL_CURRENT_TWO_MASS_CONDITIONAL_GAUSSIAN_SCALAR_CLOCK_MATTER_COUPLED_CAUSAL_INVERSE_WITH_FULL_PROFILE_ONCE_STRONG_FINITE_BALL_REMAINDER_AND_COMPLETE_CURRENT_FORCE_RECOVERY",
        "date": "2026-09-13",
        "status": "ACTUAL_CONDITIONAL_TWO_MASS_COUPLED_INVERSE_ON_SMOOTH_PREPARED_FOURIER_BALL_SUBGRAPH; NOT_UNRESTRICTED_GRAPH_SMALLNESS_STABILITY_INTERACTING_LOOPS_NONLINEAR_UV_REGGE_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/assembly.md",
            "notes/ward.md",
            "notes/remainder.md",
            "notes/inverse.md",
            "notes/recovery.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_current_assembly": serialize(
            {
                name: payload(data)
                for name, data in packets.items()
                if name
                not in (
                    "complete_actual_strong_normal_form",
                    "ordered_current_coupled_inverse",
                )
            }
        ),
        "actual_normal_form_and_inverse": serialize(
            {
                name: payload(data)
                for name, data in packets.items()
                if name
                in (
                    "complete_actual_strong_normal_form",
                    "ordered_current_coupled_inverse",
                )
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The actual QG2-H8A420 coefficient operator and complete Proca-plus-heavy conditional Gaussian response admit a two-sided scalar/clock/M1 causal inverse on each smooth prepared bounded-Fourier scalar-amplitude subgraph. The full heavy profile is matched to the current coefficient action and included exactly once; the pure Gaussian one-current, both ordered Ward legs and distinct nonlinear clock contact remain. The complete actual mode and fixed-dimensional singular-extension argument gives a strong weak-log remainder after the mixed output primitives. The reference inverse uses the entire same-prescription two-mass sum, not separate inverses or the old isolated Proca pole. An ordered weighted Volterra construction, with the variable eta-pivot commutator retained, proves smooth prepared regularity and a finite unevaluated inverse bound containing every physical kappa factor. Both original scalar constraints follow from the full weighted residual adjoint and original zero germ. The complete current force interface retains the rank-two direct auxiliary block, both time-leg clock maps and output density. The separate coherent heavy mean keeps its physical KG/source normalization. This is not an unrestricted completed-graph bound, physical smallness or stability, an interacting quantum-state theorem, nonlinear same-state bounce, quantum gravitational limit, physical UV/Regge or original V/G/B/P8 closure.",
        "not_established": [
            "A uniform bound as the external Fourier radius tends to infinity or unrestricted completed-graph surjectivity",
            "Numerically small physical inverse/forcing errors, full-gravity stability or a finite inhomogeneous nonlinear bounce",
            "Interacting light/mixed state and omitted-loop control or quantization commuting with the gravitational limit",
            "A physical cutoff, exact UV S matrix, finite-gravity Regge remainder or all-parent exclusion",
            "An unrestricted infrared shift-vector norm, literal homogeneous scalar constraints or reset initial data",
            "Original V/G/B/P8 closure or machine formalization of the continuum arguments",
        ],
        "verification_boundary": "Exact whole-profile, mixed-principal, both-leg Ward, fixed-dimensional matching, local row-primitive, ordered inverse and force identities support the written mode/singular-extension and Volterra arguments. Independent unexpanded density, full saddle, unequal-time clock/density, prepared-boundary, noncommuting matrix and weak-log diagnostics test the full formulas without replacing continuum proofs. The actual scalar all-order comparison is supported by Olbermann0704.2986v2 with the model-specific normalization and assembly proved here. Native/direct/ordinary/CLI use original SymPy; only full regression uses the separately audited exact-GCD adapter. Every frozen scientific byte remains unchanged. The continuum proofs are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The actual current two-mass coupled inverse report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.248 actual conditional two-mass coupled inverse replay passed; unrestricted physics and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
