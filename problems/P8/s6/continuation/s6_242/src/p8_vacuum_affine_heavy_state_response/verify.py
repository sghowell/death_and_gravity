"""Read-only full prepared heavy-state selection response difference and graph certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_heavy_clock_quadratic import verify as parent
from p8_vacuum_affine_heavy_curved_state import verify as state_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-heavy-state-response.json"
PARENT_SHA = "f7d3b2ddfaf2d8fb2e09dcaf28ca942ddcda8c0ac36ad9e95711b2db9a1e5a04"
STATE_SHA = "d4a6081d0482b2168255299cdd0080c90a1716eeffbc70f2645e6e2452e78eda"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_heavy_state_response/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    for prior, expected in ((parent, PARENT_SHA), (state_input, STATE_SHA)):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A frozen QG2 state, profile or full local-response input changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_241_complete_current_quadratic_and_finite_local_response_fully_rebuilt": PARENT_SHA,
        "S6_240_actual_exact_SLE_full_comparison_bounds_and_fixed_profile_fully_rebuilt": STATE_SHA,
        "all_frozen_original_and_separate_model_scientific_bytes_unchanged": True,
        "full_state_difference_not_full_comparison_response_or_quantum_inverse": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A whole state-selection response or graph gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.242.COMPLETE_EXACT_HEAVY_SLE_STATE_SELECTION_RESPONSE_DIFFERENCE_WITH_FULL_ADM_CONTACTS_ORDERED_REFLECTED_KUBO_BRIDGE_ALL_MOMENTUM_BOUNDS_AND_UNCHANGED_FIXED_PROFILE_COMMON_CLOCK_PULLBACK",
        "date": "2026-09-12",
        "status": "COMPLETE_PREPARED_STATE_SELECTION_RESPONSE_DIFFERENCE_AND_MATCHED_PROFILE; NOT_FULL_COMPARISON_STATE_RESPONSE_QUANTUM_INVERSE_INTERACTING_BOUNCE_UV_REGGE_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/hamiltonian.md",
            "notes/ordering.md",
            "notes/state.md",
            "notes/bounds.md",
            "notes/profile.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_full_ADM_covariance_and_ordered_response": serialize(
            {
                "features": payload(audit.vertices.data()),
                "response": payload(audit.response.data()),
            }
        ),
        "complete_all_momentum_state_selection_and_clock_bounds": serialize(
            payload(audit.bounds.data())
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The existing exact S240 SLE and the exact W6-initialized KG comparison define a COMPLETE prepared-state selection response DIFFERENCE, with both Cauchy covariances fixed on later histories. The full scalar ADM Hamiltonian gives every lapse, shift and noncommuting spatial-metric first/second feature. A canonical covariance tangent independently matches the full Wick commutator, including both adjoint source blocks and the instantaneous second-Hamiltonian contact. The reversed ordered term reflects both internal momenta; odd shift cross channels retain their imaginary Fourier kernels and cannot be replaced by the same ordinary imaginary part as even channels. Exact full Bogoliubov/readout inequalities, both complete pair products, full radial integrals and large-small transfer control prove all-internal/all-external-momentum bounds. The corresponding algebraic summand of the EXISTING fixed reference profile is retained without changing its sum or any prescription. Its nonzero second-clock contact cancels the current one-point chart contact only in the ENTIRE matched summand. The complete normalized state-selection-plus-profile response is below10^-1280 in the stated full physical ADM and actual common-clock weak graphs, with no external time derivative required. The full two-leg frequency projection has same-clock error below10^-1180/K on the same graph for K>=2sqrt(n). This result does not remove the full exact-comparison state's renormalized response, turn WKB evolution into an exact state, prove a same-space inverse, change the finite prescription or close the interacting clock, nonlinear bounce, quantum gravity limit, physical UV, finite-gravity Regge or original V/G/B/P8.",
        "not_established": [
            "The full exact-comparison-state renormalized heavy response or the full heavy quantum constraints/inverse",
            "A new Hadamard physical comparison state or global smoothness of auxiliary profile summands",
            "A same-space inverse, interacting linear stability, finite-amplitude remainder or nonlinear same-state bounce",
            "Control of light/mixed loops or interchange of quantization and the gravitational limit",
            "A physical cutoff, exact UV S matrix, finite-gravity Regge estimate or all-parent exclusion",
            "Original V/G/B or P8 closure",
        ],
        "verification_boundary": "Independent full Legendre and noncommuting matrix-exponential tests check the entire ADM features. A complete finite-lattice covariance and ordered Fourier calculation checks even and odd channels in nontrivially squeezed states; a separately evolved prepared metric IVP retains a nonzero final contact. Exact all-frequency analytic inequalities, not numerical fixtures, supply the full bounds and regulator tails. The fixed profile and nonlinear clock-contact chain rule are checked independently. Native/direct/ordinary/CLI retain original SymPy; only full regression uses the audited exact-GCD adapter. No frozen source, proof, test, report or prescription is changed; the written proofs are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete heavy state-selection response report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.242 complete heavy state-selection response replay passed; comparison response and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
