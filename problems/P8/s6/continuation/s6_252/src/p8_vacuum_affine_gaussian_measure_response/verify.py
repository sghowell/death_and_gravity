"""Read-only complete physical Gaussian functional and continuum Wick report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_coupled_gaussian_state import phase
from p8_vacuum_affine_coupled_gaussian_state import verify as parent

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-gaussian-measure-response.json"
)
PARENT_SHA = "e60b97e9a4634795ea788ac6a3d35618c938bad000cb57cc5a85a40bd43b1f5c"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_gaussian_measure_response/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen complete physical reference state changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_251_complete_current_free_state_and_entire_source_ancestry_fully_rebuilt": PARENT_SHA,
        "all_frozen_original_and_separate_model_scientific_bytes_unchanged": True,
        "same_fixed_H_Proca_preparations_QG2_profiles_and_M1_mean": True,
        "no_full_nonlocal_or_interacting_state_substituted": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "A full Gaussian measure, fixed-state kernel or Wick-boundary gate failed"
        )
    packets = audit.packets()
    current = phase.data()
    return {
        "schema": 1,
        "claim": "P8-S6.252.COMPLETE_CURRENT_PHYSICAL_REFERENCE_GAUSSIAN_MEASURE_AND_FIXED_STATE_CTP_COEFFICIENT_VERTICES_WITH_WRITTEN_TWO_CONE_CONTINUUM_WICK_NOISE_AND_TIME_RETARDED_EXTENSION_FAMILY",
        "date": "2026-09-13",
        "status": "EXACT_FULL_QUADRATIC_GAUSSIAN_MEASURE_AND_FINITE_CTP_WITH_WRITTEN_CONTINUUM_WICK_NOISE_AND_EXTENSION_FAMILY; NOT_NONLINEAR_COVARIANT_MEASURE_SELECTED_CURVED_COUNTERFUNCTIONAL_PHYSICAL_LOOP_MEAN_BOUND_NONLINEAR_UV_REGGE_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/measure.md",
            "notes/ctp.md",
            "notes/vertices.md",
            "notes/wick.md",
            "notes/grading.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_current_Gaussian_measure_and_fixed_state_kernels": serialize(
            {
                "full_current_coefficient_and_physical_background_map": current[
                    "entire_phase_coefficient_and_physical_background_substitution"
                ],
                "full_current_clock_functions": current[
                    "entire_current_clock_coefficients"
                ],
                "whole_measure": payload(
                    packets["whole_physical_auxiliary_Gaussian_measure"]
                ),
                "whole_fixed_state_kernels": payload(
                    packets["whole_fixed_state_quadratic_response_and_noise"]
                ),
            }
        ),
        "complete_CTP_vertices_and_continuum_Wick_extension_boundary": serialize(
            payload(
                packets["complete_CTP_coefficient_vertices_and_continuum_Wick_boundary"]
            )
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The same physical quadratic reference has a complete finite-regulator auxiliary Dirac/configuration Gaussian measure, retaining its nonzero lower constraint bracket, density factors and crossing-safe canonical interpretation. The fixed-state CTP trace has its full continued metaplectic phase and every coefficient first/second vertex; independent full-matrix Wick, Duhamel and overlap routes agree on means, time-retarded susceptibility, seagulls and positive connected noise. The written two-cone graph/smoothing proof supplies continuum reference Wick bilinears and their full noise. Complete chart symbol orders give a conservative finite diagonal scaling bound and a family of time-retarded distributional extensions with local ambiguities. Existing counterprofile grades and nonstationary embedding contacts remain explicit. This is not a complete nonlinear covariant/BRST measure, physical parent vertex correspondence, selected curved counterfunctional, renormalized metric/light mean, evaluated omitted-loop bound, interacting state, compatible nonlinear bounce, physical UV/Regge or original V/G/B/P8 closure.",
        "not_established": [
            "A full nonlinear covariant/BRST, ghost, connection or field-map quantum measure",
            "The complete off-reference physical metric/light vertex and variational/Ward-consistent curved counterfunctional",
            "A selected physical finite subtraction, renormalized physical metric/light mean or quantitative omitted-loop bound",
            "The full nonlocal/interacting state, convergence of perturbation theory, compatible nonlinear norm or same-state nonlinear bounce",
            "Spacelike microcausality from a gauge-reduced time-retarded kernel, physical UV scattering or finite-gravity IR/Regge control",
            "Original V/G/B/P8 closure or formalization of the continuum argument",
        ],
        "verification_boundary": "All full quadratic identities use exact arithmetic and every frozen source remains unchanged. Independent Gaussian, oscillator, coefficient and phase diagnostics support but do not FORMALIZE the written continuum Wick/scaling/extension arguments. A reference Wick mean is not a covariant physical stress. Native/direct/ordinary/CLI use original SymPy; only full regression uses the audited exact-GCD adapter.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete Gaussian functional report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.252 full reference Gaussian functional replay passed; physical curved subtraction and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
