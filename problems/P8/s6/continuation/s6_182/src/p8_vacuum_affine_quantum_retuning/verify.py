"""Read-only bounded fixed retuning of the regular affine parent."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_kernel_norm import verify as parent
from p8_vacuum_affine_proca_gaussian import verify as gaussian

from . import audit, bounds, profile, response, ward

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-quantum-retuning.json"
PARENT_SHA = "7ec845d86ab05c33addb35a51380b874b73ee73e58702bf195136c2751acb8ea"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_quantum_retuning/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA or sha(gaussian.REPORT) != GAUSSIAN_SHA:
        raise ValueError(
            "The source-pinned regular affine Gaussian/kernel chain changed"
        )
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_181_fully_rebuilt": PARENT_SHA,
        "S6_176_same_Gaussian_stress_state_scheme_transitively_rebuilt": GAUSSIAN_SHA,
        "new_action_declared_separately_not_written_over_old_parent": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A separately named fixed-retuning gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.182.NEW_REGULAR_AFFINE_QG1_FIXED_COEFFICIENT_STATIONARY_GAUSSIAN_CD_AND_VACUUM_WITH_EXPLICIT_ERROR_BUDGET_AND_ACTUAL_RESPONSE",
        "date": "2026-09-11",
        "status": "NEW_ACTION_RETAINED_GAUSSIAN_CD_AND_MINKOWSKI_MEAN_STATIONARITY_WITH_BOUNDED_COEFFICIENT_CORRECTION_AND_ACTUAL_PREPARED_RESPONSE; NOT_FULL_INTERACTING_QUANTUM_BACKGROUND_STABILITY_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/profile.md",
            "notes/stationary.md",
            "notes/bounds.md",
            "notes/response.md",
            "notes/vacuum.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "new_fixed_coefficient_and_stationary_mean": serialize(
            {"profile": payload(profile.data()), "ward": payload(ward.data())}
        ),
        "quantitative_error_budget_and_actual_response": serialize(
            {"bounds": payload(bounds.data()), "response": payload(response.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The separately named CD-REG-AFFINE-ISO-QG1 adds an explicit fixed scalar coefficient, computed once from the specified reference Gaussian stress, without changing the state, vector operator or covariant scheme. It cancels the retained CD one-point currents and scalar Ward source and the same-prescription Minkowski vacuum energy. All36 mixed coefficient derivatives through five in each variable are bounded across the full original X strip and more sharply near the clock. Nonconstant vacuum field jets below2048 vanish, and the complete anchored canonical interactions are fixed. The actual stationary response retains the nonlocal vector scale channel and its matched prepared inverse, with local lapse C1 below15000. This is not a full interacting quantum background, small/stable full inverse, physical cutoff or original V/G/B closure.",
        "not_established": [
            "A full interacting affine/light/tensor/auxiliary/mixed-loop state or quantum background",
            "A useful small full inverse, nonlinear stability or spatial/noise/global-time perturbation control",
            "Physical heavy/cutoff matching, omitted-order and threshold errors across the common parent domain",
            "Complete vacuum scattering cuts/contour/truncation, finite-gravity IR/Regge remainder or original V/G/B closure",
        ],
        "verification_boundary": "Exact complete coefficient/jet/covariant normalization/Ward/current/chart/canonical-family identities, continuous complex-disc Cauchy estimates and actual local constraint/response bounds, with independent literal and high-precision fixtures. The reference stress is fixed by the rebuilt all-order source chain, not freely assumed. Native, direct science, ordinary and CLI use original SymPy; only full regression uses its audited exact GCD adapter. Written analytic and Hadamard arguments are not FORMALIZED. A new finite coefficient choice is not a UV-matching certificate.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The new fixed quantum-retuning report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.182.NEW_REGULAR_AFFINE_QG1_FIXED_COEFFICIENT_STATIONARY_GAUSSIAN_CD_AND_VACUUM_WITH_EXPLICIT_ERROR_BUDGET_AND_ACTUAL_RESPONSE replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
