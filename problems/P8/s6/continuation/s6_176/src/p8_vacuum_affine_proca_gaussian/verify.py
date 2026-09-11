"""Read-only replay of the actual conditional Gaussian vector state and stress."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_proca_stress import verify as stress_parent
from p8_vacuum_affine_vector_schur import verify as parent
from p8_vector_hadamard import verify as state_parent

from . import audit, bridge, gaussian, state, stress

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-proca-gaussian.json"
PARENT_SHA = "151f0642185a494dd1a49e77360349c5105596fbf4485cdb911586ef0d480b94"
STATE_SHA = "8ba1aa16cee4dbe890b233fac864a83de569402e06b26c9c18ae56ef763a9200"
STRESS_SHA = "4e3c6103906b8ddb344806b7fe9309493f3c8dfda28629e70912872be65c09f8"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_proca_gaussian/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    for module, digest in (
        (parent, PARENT_SHA),
        (state_parent, STATE_SHA),
        (stress_parent, STRESS_SHA),
    ):
        if sha(module.REPORT) != digest:
            raise ValueError(
                "A frozen Gaussian-parent or explicit state/stress source changed"
            )
        module.validate_report(
            json.loads(module.REPORT.read_text()), module.build_report()
        )
    return {
        "S6_175_fully_rebuilt": PARENT_SHA,
        "S6_55_explicit_all_order_state_rebuilt": STATE_SHA,
        "S6_82_ordinary_Proca_stress_rebuilt": STRESS_SHA,
        "same_classical_target_and_all_prior_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A conditional Gaussian vector proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.176.SOURCE_AWARE_CONDITIONAL_AFFINE_PROCA_HADAMARD_FAMILY_AND_FIXED_CLOCK_GAUSSIAN_STRESS",
        "date": "2026-09-11",
        "status": "CONDITIONAL_SOURCED_VECTOR_GAUSSIAN_HADAMARD_FAMILY_AND_SPECIFIED_CLOCK_C5_STRESS; NOT_FULL_PARENT_QUANTUM_SOLUTION_MATCHING_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/bridge.md",
            "notes/state.md",
            "notes/renormalization.md",
            "notes/stress.md",
            "notes/response.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "conditional_Gaussian_state_and_canonical_bridge": serialize(
            {
                "bridge": payload(bridge.data()),
                "clock": payload(bridge.clock()),
                "state": payload(state.data()),
            }
        ),
        "specified_vector_stress_and_reference_clock_response": serialize(
            {
                "stress": payload(stress.data()),
                "prescription": payload(stress.prescription()),
                "local": payload(stress.local()),
                "Gaussian": payload(gaussian.data()),
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The same regular classical affine parent admits an explicitly normalized conditional Gaussian vector sector at zeta=1e-6 and kappa=1e800. Its actual arbitrary-lapse/scale readouts and complete clock operator equal ordinary mass1000 Proca, with a stated canonical measure, fixed all-order Cauchy covariance and covariant finite prescription. Compact sourced histories have a positive coherent Hadamard family with unchanged connected Proca fluctuations. On the reference clock the twelve specified energy/pressure time derivatives through order five are below1e30 on [-1/2,1/2], or1e-770 against kappa; no old scalar stress-canceling profiles are installed. The source-induced mean force starts cubically. These conditional-sector identities and bounds do not supply a full interacting parent, a self-consistent quantum bounce, physical cutoff, full real-time matching or V/G/B closure.",
        "not_established": [
            "Quantitative nonlinear sourced mean, symmetric noise or functional metric-response norms",
            "Full affine/metric/light/tensor quantum measure, auxiliary and mixed loops or a controlled interacting state",
            "A self-consistent quantum-corrected background, physical cutoff or all omitted-order bounds",
            "Common vacuum matching, contour/truncation, finite-gravity IR/Regge remainder, V/G/B or original P8 closure",
        ],
        "verification_boundary": "Exact canonical, physical metric variation, covariance, source-degree, local Ward and normalization identities bridge fully replayed source-pinned ordinary-Proca calculations. Independent lapse/scale variations, constrained canonical evolution, state/normalization and source-contact controls supplement written coherent-state and theorem-hypothesis proofs. This is not formalization or a full quantum-parent certificate. Native, direct science, ordinary and CLI use original SymPy; full regression alone uses the separately audited exact GCD adapter.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The conditional affine-Proca Gaussian report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.176.SOURCE_AWARE_CONDITIONAL_AFFINE_PROCA_HADAMARD_FAMILY_AND_FIXED_CLOCK_GAUSSIAN_STRESS replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
