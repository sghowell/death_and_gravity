"""Read-only actual-CD source-force covariance and canonical-normalization report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_proca_gaussian import verify as gaussian
from p8_vacuum_affine_retarded_energy import verify as retarded
from p8_vacuum_affine_vacuum_causal_noise import verify as parent

from . import audit, covariance, derivatives, force, modes

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-cd-source-noise.json"
PARENT_SHA = "2800f43a4e682b4751bffccd91d999ff6d3e54994b186fdf7c9009cea027a91b"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
RETARDED_SHA = "c25ba2599f6297513798ff45291a831aa46be8397dad43fec3d6dc8fdfaea488"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_cd_source_noise/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    for previous, expected in (
        (parent, PARENT_SHA),
        (gaussian, GAUSSIAN_SHA),
        (retarded, RETARDED_SHA),
    ):
        if sha(previous.REPORT) != expected:
            raise ValueError(
                "The source-pinned affine state and full source chain changed"
            )
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_184_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "S6_178_complete_clock_source_transitively_rebuilt": RETARDED_SHA,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("An actual CD source-noise gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.185.ACTUAL_ALL_ORDER_CD_STATE_COMPLETE_AFFINE_SOURCE_SMEARED_SCALAR_FORCE_COVARIANCE_AND_EXPLICIT_CANONICAL_NORMALIZATION",
        "date": "2026-09-11",
        "status": "SAME_ALL_ORDER_CD_STATE_FULL_SOURCE_SCALAR_FORCE_NOISE_BOUND_WITH_CANONICAL_KAPPA_CANCELLATION; NOT_METRIC_NOISE_QUANTUM_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/state.md",
            "notes/covariance.md",
            "notes/source-gradient.md",
            "notes/force.md",
            "notes/normalization.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "actual_state_and_curved_covariance": serialize(
            {"modes": payload(modes.data()), "covariance": payload(covariance.data())}
        ),
        "full_source_derivatives_and_force_conventions": serialize(
            {"derivatives": payload(derivatives.data()), "force": payload(force.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The actual unchanged all-order CD Gaussian state admits a continuous full three-polarization covariance upper bound, including the physical temporal constraint, evolving canonical normalizers and spatial gradients. A full complex-jet proof controls the complete near-clock source's mixed spatial Frechet derivative. The clock-normalized force has standard deviation below2e-392 delta per unit raw-clock H3 test norm, but canonical normalization cancels this gravitational suppression and gives the bound2e8 delta. On the explicit smaller ball delta<=1e-14 the canonical source-noise bound is below2e-6 and the previous full-source mean relative coefficient is below4e-18. Reference source-induced noise vanishes exactly; metric stress noise does not follow. This is a scoped conditional scalar observable, not a quantum solution, physical cutoff or original V/G/B closure.",
        "not_established": [
            "Metric stress-tensor covariance, coincident noise, spatial/nonlinear/stochastic light or gravity solutions",
            "Small or stable full response, finite-coupling truncation error or an interacting quantum background",
            "Interacting parent measures/loops, complete heavy/cutoff/threshold and omitted-order matching",
            "Complete vacuum cuts/contour/truncation, finite-gravity IR/Regge or original V/G/B closure",
        ],
        "verification_boundary": "Exact actual mode/normalizer/constraint and full source identities; source-pinned all-order preparation and evolution envelopes; continuous complex Cauchy, all-momentum covariance and Sobolev bounds. Independent actual all-order Cauchy, time-evolved mode, complex full-switch, covariant-source and canonical-normalization fixtures supplement the written proofs. Native, direct science, ordinary and CLI use original SymPy; only full regression uses the audited exact GCD adapter. Analytic continuum arguments are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The actual CD source-noise report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.185.ACTUAL_ALL_ORDER_CD_STATE_COMPLETE_AFFINE_SOURCE_SMEARED_SCALAR_FORCE_COVARIANCE_AND_EXPLICIT_CANONICAL_NORMALIZATION replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
