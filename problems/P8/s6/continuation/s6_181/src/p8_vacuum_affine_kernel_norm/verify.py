"""Read-only numeric full-half-line bound on the same massive scalar kernel."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_proca_rank_one_inverse import verify as old_scalar
from p8_vacuum_affine_coupled_response import verify as parent

from . import audit, dyadic, estimate, tail, threshold

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-kernel-norm.json"
PARENT_SHA = "ad283ebc8d7159256c851b51a1f46056ec5ccb1924e34e6f6f02a28b0866aa3d"
SCALAR_SHA = "0f04fd277040ad21525eb1d9102e6af45c4d695fac3f69d0d37710a7ab8a0db6"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_kernel_norm/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA or sha(old_scalar.REPORT) != SCALAR_SHA:
        raise ValueError("The frozen current inverse or scalar density source changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_180_fully_rebuilt": PARENT_SHA,
        "S6_86_same_scalar_kernel_transitively_rebuilt": SCALAR_SHA,
        "new_continuous_constants_derived_from_literal_source_density": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("An explicit massive scalar-kernel norm gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.181.SAME_MASSIVE_SCALAR_RANGE_KERNEL_EXPLICIT_MASS_INDEPENDENT_FULL_HALF_LINE_L1_BOUND",
        "date": "2026-09-11",
        "status": "EXPLICIT_ALL_FREQUENCY_SCALAR_KERNEL_NORM_BELOW_3E11_FOR_EVERY_POSITIVE_MASS; NOT_SMALL_FULL_INVERSE_QUANTUM_BACKGROUND_STABILITY_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/threshold.md",
            "notes/tail.md",
            "notes/partition.md",
            "notes/norm.md",
            "notes/transfer.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "continuous_density_and_band_bounds": serialize(
            {
                "threshold": payload(threshold.data()),
                "tail": payload(tail.data()),
                "dyadic": payload(dyadic.data()),
            }
        ),
        "same_kernel_norm_and_current_inverse_transfer": serialize(
            payload(estimate.data())
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The same normalized massive scalar range inverse has full half-line L1 norm below201600000000 and hence below3e11 for every positive mass. Explicit continuous endpoint/tail derivative envelopes and a C2 dyadic partition control every frequency band. Absolute L1 summation and the nonnegative spectral Laplace representation identify the same inverse, without discarded poles, instantaneous terms or tails. This replaces only the K_L1 parameter in the current retained coupled inverse; the full curved weak-log majorant C remains unevaluated. It does not yield a small or stable full inverse, a quantum background or original P8 closure.",
        "not_established": [
            "A numeric full curved weak-log remainder majorant or small/stable coupled inverse",
            "A compatible quantum preparation and controlled nonlinear corrected background, spatial/noise or global-time response",
            "Full interacting measure/state/loops, physical cutoff or omitted-order error control",
            "Complete vacuum matching/contour/cuts/truncation, finite-gravity IR/Regge or original V/G/B closure",
        ],
        "verification_boundary": "Exact density/derivative/product/partition/normalization identities and continuous written endpoint/tail/Fourier/Laplace estimates, with independent high-precision density, derivative and spectral-transform checks. The source-pinned current inverse chain is rebuilt unchanged. Native, direct science, ordinary and CLI use original SymPy; only full regression uses its audited exact GCD adapter. Not FORMALIZED, an optimal norm or a full quantum-parent certificate.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The scalar-kernel norm report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.181.SAME_MASSIVE_SCALAR_RANGE_KERNEL_EXPLICIT_MASS_INDEPENDENT_FULL_HALF_LINE_L1_BOUND replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
