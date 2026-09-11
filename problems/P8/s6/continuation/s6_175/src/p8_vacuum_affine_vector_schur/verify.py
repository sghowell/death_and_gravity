"""Read-only replay of exact sourced-vector reduction and scoped estimates."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_analytic_affine_parent import verify as parent

from . import audit, gaussian, homogeneous, lorentz, vacuum

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-vector-schur.json"
PARENT_SHA = "e0cf546e77cf12ff4e0f2e6bef706ff6de1eddc40cacff3b0b151fcb7009a4f4"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_vector_schur/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen regular classical affine-parent report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_174_fully_rebuilt": PARENT_SHA,
        "same_regular_classical_affine_parent_and_complete_target": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A sourced-vector proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.175.EXACT_SOURCED_AFFINE_VECTOR_SCHUR_REDUCTION_AND_SCOPED_VACUUM_RESPONSE_BOUNDS",
        "date": "2026-09-11",
        "status": "EXACT_CLASSICAL_LIFT_AND_COMMON_REGULATOR_SCHUR_FULL_EUCLIDEAN_AND_LEADING_LORENTZ_CONTROL; NOT_FULL_SOURCE_REAL_TIME_QUANTUM_MATCHING_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/operator.md",
            "notes/homogeneous.md",
            "notes/gaussian.md",
            "notes/vacuum.md",
            "notes/lorentz.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "exact_sourced_vector_and_homogeneous_reduction": serialize(
            {
                "operator": payload(gaussian.data()),
                "variation": payload(gaussian.variation()),
                "homogeneous": payload(homogeneous.data()),
            }
        ),
        "full_source_and_finite_germ_estimates": serialize(
            {"vacuum": payload(vacuum.data()), "Lorentz": payload(lorentz.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "For the same separately named CD-REG-AFFINE-ISO parent, the full source gives an exact closed-source/homogeneous classical lift, a causal fixed-background sourced Proca inverse and a common-regulator Gaussian Schur identity with all contact and divergence terms retained. The full analytic source has an explicit small-field L2 bound and a tiny independent Euclidean source-action bound. Its leading vacuum degree-four source has a separate below-pole Lorentzian estimate using the complete propagator. The full nonlinear source is not assumed bandlimited; no complete real-time/state/influence bound, renormalized interacting stress, physical cutoff, V/G/B or original P8 closure follows.",
        "not_established": [
            "A quantitative full-source real-time matching or in-in influence bound on the stated vacuum class",
            "A selected vector Hadamard state, renormalized stress/metric response or controlled interacting quantum background",
            "A physical cutoff, all omitted-order bounds or transfer of old GY14 matching to this parent",
            "Vacuum contour/truncation, finite-gravity IR/Regge remainder, V/G/B or original P8 closure",
        ],
        "verification_boundary": "Exact matrix/operator identities, explicit full-function inequalities and written causal/functional proofs are supplemented by independent coordinate Euler equations, dense Gaussian solves, high-precision actual-source controls, support/pole and source-omission tests. This is not a formalized theorem or a renormalized quantum-state certificate. Native, direct science, ordinary and CLI use original SymPy; full regression alone uses its separately audited exact GCD adapter.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The sourced affine-vector Schur report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.175.EXACT_SOURCED_AFFINE_VECTOR_SCHUR_REDUCTION_AND_SCOPED_VACUUM_RESPONSE_BOUNDS replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
