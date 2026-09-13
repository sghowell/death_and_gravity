"""Read-only full physical Gaussian reference and reduced short-distance report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_heavy_clock_quadratic import verify as clock_input
from p8_vacuum_affine_heavy_source_filtration import verify as parent
from p8_vacuum_affine_two_mass_nonscalar_inverse import verify as tensor_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-coupled-gaussian-state.json"
PARENT_SHA = "bc33180dd5d312f638241d8f04e4ddf6eb2863d6fcf85015c44a8ba90bcfd77c"
CLOCK_SHA = "f7d3b2ddfaf2d8fb2e09dcaf28ca942ddcda8c0ac36ad9e95711b2db9a1e5a04"
TENSOR_SHA = "619218777bc9d5e525a04d7276af6fa9fcf1f46826d7badcf35243859ba48b22"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_coupled_gaussian_state/*.py"))
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
        (clock_input, CLOCK_SHA),
        (tensor_input, TENSOR_SHA),
    ):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A frozen whole current phase, tensor or source-boundary input changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_250_source_independent_quantum_frontier_and_complete_ancestry_fully_rebuilt": PARENT_SHA,
        "S6_241_complete_current_coefficient_functions_and_tame_chart_hypotheses_fully_rebuilt": CLOCK_SHA,
        "S6_249_actual_classical_nonscalar_action_and_global_TT_projector_fully_rebuilt": TENSOR_SHA,
        "all_frozen_original_and_separate_model_scientific_bytes_unchanged": True,
        "full_quantum_response_not_replaced_by_free_state_construction": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError(
            "A full physical Gaussian reference or short-distance input gate failed"
        )
    packets = audit.packets()
    return {
        "schema": 1,
        "claim": "P8-S6.251.COMPLETE_CURRENT_PHYSICAL_COUPLED_SCALAR_AND_TENSOR_FIXED_GAUSSIAN_REFERENCE_WITH_FULL_CANONICAL_CCR_MATRIX_SAMPLING_AND_WRITTEN_ALL_ORDER_REDUCED_TWO_CONE_SHORT_DISTANCE_CONTROL",
        "date": "2026-09-13",
        "status": "EXACT_LOCAL_PHYSICAL_QUADRATIC_REFERENCE_GAUSSIAN_STATE_WITH_WRITTEN_REDUCED_TWO_CONE_SHORT_DISTANCE_PROOF; NOT_FULL_NONLOCAL_INTERACTING_STATE_COVARIANT_COUNTERFUNCTIONAL_LOOP_SIZE_NONLINEAR_UV_REGGE_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/state.md",
            "notes/ultraviolet.md",
            "notes/constraints.md",
            "notes/tensor.md",
            "notes/subtraction.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_physical_phase_and_preparation": serialize(
            payload(packets["complete_current_canonical_phase_and_preparation"])
        ),
        "complete_Gaussian_state_and_short_distance_input": serialize(
            {
                name: payload(data)
                for name, data in packets.items()
                if name != "complete_current_canonical_phase_and_preparation"
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The full current local coefficient Hamiltonian has a newly specified physical free reference for its two coupled scalar and two tensor modes. Canonical variables retain both kappa and volume factors, complete Legendre maps, gyroscopic terms and time-dependent chart boundaries. A fixed smooth positive preparation Gramian defines a basis-independent pure Gaussian covariance with exact complex positivity, CCR and unequal-time propagation; its full mixing is retained and the zero-momentum extension is infrared finite. The tensor covariance uses the same preparation rule and bounded coordinate-free TT projector. The written all-order homogeneous matrix proof retains both current characteristic cones, every subleading symbol coefficient, positive-sum Riccati recursion, exact polynomial-loss comparison and full anomalous-block smooth sampling across the chart cover. It supplies the reduced oriented short-distance condition, not a scalar-theorem shortcut or finite-order formalization. The state is fixed under later variation. Existing Proca/H states and the classical M1 mean remain unchanged. This local physical free reference is not the full nonlocal/interacting quantum state, a covariant gauge/constraint measure, a curved counterfunctional or a renormalized physical stress/omitted-loop estimate. No new mean retuning, compatible nonlinear bounce, quantum gravitational limit, physical UV/Regge or original V/G/B/P8 closure follows.",
        "not_established": [
            "A full nonlocal conditional or interacting light/metric/M1 quantum state",
            "Covariant gauge/constraint or ghost measure, full curved counterfunctional or renormalized physical metric/light stress",
            "Numerical omitted-loop bounds, convergence of perturbation theory or a complete quantum gravitational limit",
            "A nonlinear same-state bounce, compatible evaluated feedback norm, physical smallness or stability",
            "Global-time completeness of the interacting theory, all-parent exclusion, physical UV scattering or finite-gravity IR/Regge estimate",
            "Original V/G/B/P8 closure or machine formalization of the continuum proof",
        ],
        "verification_boundary": "Full current Hamiltonian, weighted/canonical CCR, complete chart and boundary transformations, noncommuting Gaussian covariance, complex uncertainty, actual principal pencils and full matrix Riccati identities support the written reduced free-reference state proof. Independent exact and numerical phase, symplectic covariance, integrated preparation and complete-mode diagnostics do not replace the arbitrary-order symbol argument or construct the missing interacting theory. Native/direct/ordinary/CLI use original SymPy; only full regression uses the audited exact-GCD adapter. Frozen scientific bytes remain unchanged. These arguments are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The whole physical Gaussian reference report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.251 full physical free-reference state replay passed; interacting physics and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
