"""Read-only new QG2 clock comparison and complete finite heavy local response report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_heavy_curved_state import verify as parent
from p8_vacuum_affine_heavy_scalar_parent import verify as classical_input
from p8_vacuum_affine_scalar_tame_propagator import verify as comparison_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-heavy-clock-quadratic.json"
PARENT_SHA = "d4a6081d0482b2168255299cdd0080c90a1716eeffbc70f2645e6e2452e78eda"
COMPARISON_SHA = "7560a2ba6fcb470f784662046e2561d641bf4f6f0034104a158d977093075255"
CLASSICAL_SHA = "934b6c035fb5de0eb89ee1e1a03c2624c78a33985c0503bce883c41240b05905"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_heavy_clock_quadratic/*.py"))
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
        (comparison_input, COMPARISON_SHA),
        (classical_input, CLASSICAL_SHA),
    ):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A frozen QG2 state, full chart proof or covariant parent changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_240_complete_QG2_state_stress_and_fixed_profile_fully_rebuilt": PARENT_SHA,
        "S6_221_generic_full_classical_chart_energy_and_all_matrix_enclosures_fully_rebuilt": COMPARISON_SHA,
        "S6_238_complete_covariant_R_and_source_clock_jets_fully_rebuilt": CLASSICAL_SHA,
        "all_frozen_original_and_separate_model_scientific_bytes_unchanged": True,
        "current_quadratic_and_local_result_not_full_quantum_inverse_or_bounce": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "A new clock hypothesis, whole local Hessian or graph gate failed"
        )
    return {
        "schema": 1,
        "claim": "P8-S6.241.FULL_QG2_CLOCK_QUADRATIC_WITH_EXPLICIT_CURRENT_COMPARISON_HYPOTHESES_AND_COMPLETE_PHYSICAL_HEAVY_FINITE_LOCAL_RESPONSE_AND_FIXED_PROFILE_GRAPH_BOUNDS",
        "date": "2026-09-12",
        "status": "CURRENT_CLASSICAL_COEFFICIENT_COMPARISON_AND_COMPLETE_FINITE_LOCAL_HEAVY_RESPONSE; NOT_FULL_QUANTUM_RESPONSE_INVERSE_INTERACTING_BOUNCE_UV_REGGE_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/clock.md",
            "notes/propagator.md",
            "notes/local.md",
            "notes/tensor.md",
            "notes/bounds.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_current_clock_quadratic_comparison": serialize(
            payload(audit.clock.data())
        ),
        "complete_physical_scalar_and_tensor_local_response": serialize(
            {
                "scalar": payload(audit.local.data()),
                "tensor": payload(audit.tensor.data()),
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The full existing QG2 profile and covariant R are expanded in the physical-to-affine lapse chart, retaining the pivot, mixed and volume terms. The actual combined Proca-plus-H fixed reference stress is explicitly shown to satisfy every generic S221 chart/coefficient hypothesis; no old stability conclusion is transferred automatically. The full current classical coefficient-sector phase comparison retains exp(1e29)(1+P²)^6 with12 spatial derivatives lost. The independent added classical heavy mode has a full exact energy-scaled propagation bound below1000. Its Gaussian metric response is not zero. The ENTIRE finite scalar heat action and fixed heavy-profile quadratic are computed in physical ADM lapse, metric and shift variables and pulled back to the actual common clock, with complete second-metric, curvature and nonzero state-profile chart-contact terms. The constant heavy vacuum action cancels covariantly before any variation or constraint elimination. The matched local heat block plus the retained fixed state-profile Hessian has complete physical AND same-clock normalized weak-graph bounds below10^-580 on[-1/2,1/2], uniformly over all momentum. The full tensor finite Euler operator retains its actual a³ pairing, variable A_H and every fourth-order contact, and obeys a10^-597 bound on the specified order-four graph. These are complete LOCAL response and classical coefficient-comparison results, not the remaining exact SLE/subtraction-dependent determinant response or its full coupled inverse. Interacting light/mixed clock loops, quantum gravitational decoupling, same-state nonlinear bounce, physical UV, finite-gravity Regge and original V/G/B/P8 remain open.",
        "not_established": [
            "The remaining full state/subtraction-dependent heavy determinant response or a complete quantum lapse/shift inverse",
            "A same-space classical or quantum inverse from coefficient smallness, or full interacting linear or nonlinear stability",
            "An interacting light/metric clock state, control of omitted loops or a quantum gravitational-limit theorem",
            "Global-time numerical bounds beyond the named intervals, an interpolating history or a same-state nonlinear bounce",
            "A physical interpretation of isolated local-operator poles, a physical cutoff, exact UV S matrix or finite-gravity Regge bound",
            "All-parent exclusion or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact full metric/curvature/profile variations are independently checked against a nonlinear ADM identity, whole-action lapse/scale Euler equations and a flat Weyl/Euler-boundary contraction. The physical-to-clock pullback checks the full second-variation chain rule: the matched heat block is on shell, but the state-profile contact alone is nonzero and retained. Exact current coefficient envelopes license the generic full two-chart proof; its source-pinned algebra is not represented as a new theorem. Complete local coefficient sums give the physical AND same-clock graphs, with every shift, volume, second-metric and weighted-adjoint term retained. Numerical tests are diagnostics, not full quantum-response estimates. Native/direct/ordinary/CLI retain original SymPy; only full regression uses the audited exact-GCD adapter. No frozen scientific or prescription bytes change; these written proofs are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete new clock-comparison and finite local-response report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.241 current clock comparison and complete finite local response replay passed; full quantum inverse and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
