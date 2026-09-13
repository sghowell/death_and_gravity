"""Read-only specified heavy Gaussian state and retained quantum mean report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_heavy_parent_one_loop import verify as parent
from p8_vacuum_affine_proca_gaussian import verify as gaussian_input
from p8_vacuum_affine_quantum_retuning import verify as original_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-heavy-curved-state.json"
PARENT_SHA = "a09852126d79ca223af06555e10283328cb144a14945c127abee73b7ce8b9c05"
ORIGINAL_SHA = "f754f6d62795d5a1c0ade4c576395f1587fc35f9c9895bc8b7b15b2b2dec2e16"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_heavy_curved_state/*.py"))
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
        (original_input, ORIGINAL_SHA),
        (gaussian_input, GAUSSIAN_SHA),
    ):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A frozen full parent, vacuum loop or original Gaussian input changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_239_complete_formal_classical_limit_first_loop_fully_rebuilt": PARENT_SHA,
        "S6_182_original_fixed_QG1_reference_profile_and_vacuum_matching_fully_rebuilt": ORIGINAL_SHA,
        "S6_176_retained_original_Proca_state_and_prescription_fully_rebuilt": GAUSSIAN_SHA,
        "all_frozen_original_and_separate_model_scientific_bytes_unchanged": True,
        "new_heavy_state_and_profile_not_full_interacting_quantum_parent": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "A full heavy-state, dimensional, stress or mean-profile gate failed"
        )
    return {
        "schema": 1,
        "claim": "P8-S6.240.SPECIFIED_HEAVY_SLE_HADAMARD_STATE_WITH_FULL_COVARIANT_SCALAR_RENORMALIZATION_COMPLETE_ALL_MOMENTUM_REFERENCE_STRESS_BOUNDS_AND_FIXED_QG2_MEAN_PROFILE_PRESERVING_FIRST_LOOP_VACUUM_MATCHING",
        "date": "2026-09-12",
        "status": "SPECIFIED_GAUSSIAN_HEAVY_STATE_REFERENCE_STRESS_AND_RETAINED_MEAN_MATCHING; NOT_FULL_INTERACTING_CLOCK_RESPONSE_NONLINEAR_BOUNCE_UV_REGGE_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/state.md",
            "notes/renormalization.md",
            "notes/bounds.md",
            "notes/profile.md",
            "notes/scope.md",
            "notes/literature.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "actual_heavy_state_full_covariant_stress_and_fixed_mean_profile": serialize(
            {
                "state": payload(audit.state.data()),
                "quantum": payload(audit.quantum.data()),
            }
        ),
        "complete_uniform_stress_and_new_profile_bounds": serialize(
            payload(audit.estimates.data())
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "A fixed compact sampling before t0=-1/2 specifies the added heavy scalar's exact state of low energy on the full CD reference. Physical positive frequency, CCR, the regular energy minimizer and basis independence are fixed; the applicable primary theorem supplies Hadamard regularity. A separate explicit analytic estimate retains the whole WKB comparison residual, full state-selection error, all momentum and the complete fourth-order subtraction. The dimensionally continued scalar pole action is varied before the limit, yielding the full finite mu1 matching rather than component-only subtraction. Both reference stress components and all time derivatives0..5 are below10^-400 against kappa0 on[-1,1]. A separately named QG2 physical scalar coefficient, fixed once from the actual source-defined stress, restores the retained Proca-plus-H Gaussian mean equations on the same physical CD metric with the recomputed S238 affine mean. Heavy and light flat Gaussian vacuum constants are explicitly retained and cancelled in the formal limiting-action one-loop vacuum. The nonconstant profile starts at field degree2048, so the S239 complete first-loop four-Phi result remains unchanged. The combined full normalized clock four-jets fit a NEW10^-399 budget on the named tube. The new scalar clock jets change the light Hessian: old stability and quantum inverse results do not automatically transfer. Full interacting light/mixed clock loops and state, quantum gravitational decoupling, full new response/inverse, nonlinear same-state bounce, physical UV, finite-gravity Regge and original V/G/B/P8 remain open.",
        "not_established": [
            "A full interacting light/metric quantum state on the clock or control of all light/mixed and omitted loops",
            "A full new quantum response or compatible inverse bound, automatic old linear-stability transfer, or nonlinear same-state bounce",
            "A quantitative global-time bound outside[-1,1], or global real analyticity of the smooth reference profile",
            "Quantum gravitational decoupling, a quantum-vacuum stability theorem, physical UV completion or finite-gravity Regge estimates",
            "A physical cutoff, all-parent exclusion or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact physical mode, CCR, dimensional-counteraction and profile identities accompany written complex-domain, full Volterra, repeated integration-by-parts and absolute radial-integral proofs with explicit constants. Independent finite-mass numerical exact-mode fixtures test the state construction, not the actual enormous hierarchy. The Hadamard theorem's hypotheses are checked against the primary paper; these arguments are not FORMALIZED. Native/direct/ordinary/CLI retain original SymPy; only full regression uses the audited exact-GCD adapter. No frozen scientific or prescription bytes are changed.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete specified heavy-state and reference-mean report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.240 specified heavy Gaussian state and retained mean replay passed; full interacting clock and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
