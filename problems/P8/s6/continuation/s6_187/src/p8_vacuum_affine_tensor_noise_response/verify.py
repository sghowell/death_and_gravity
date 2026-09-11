"""Read-only restricted, tree-propagated CD tensor-noise response bound."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_cd_metric_noise import verify as parent
from p8_vacuum_affine_proca_gaussian import verify as gaussian

from . import audit, detectors, energy, projector, tree

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-tensor-noise-response.json"
PARENT_SHA = "cdbcd5e9a070b1c1c4d24ebc4f0b602ef71d8ddce6f85381cf7cb55e44343f2c"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_tensor_noise_response/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    for previous, expected in ((parent, PARENT_SHA), (gaussian, GAUSSIAN_SHA)):
        if sha(previous.REPORT) != expected:
            raise ValueError(
                "The source-pinned actual affine Gaussian state chain changed"
            )
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_186_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("An actual compatible tensor-response gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.187.ACTUAL_CD_BOUNDARY_COMPATIBLE_TREE_TENSOR_PROPAGATED_STRESS_NOISE_AND_ALL_MOMENTUM_ADJOINT_BOUND",
        "date": "2026-09-11",
        "status": "ACTUAL_CD_RESTRICTED_COMPATIBLE_TT_DETECTOR_TREE_PROPAGATED_NOISE_BOUND; NOT_FULL_QUANTUM_TENSOR_STATE_SELF_ENERGY_REMAINDER_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/tree.md",
            "notes/projector.md",
            "notes/boundary.md",
            "notes/energy.md",
            "notes/noise.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "actual_tensor_action_and_TT_class": serialize(
            {"tree": payload(tree.data()), "projector": payload(projector.data())}
        ),
        "adjoint_energy_and_compatible_noise": serialize(
            {"energy": payload(energy.data()), "detectors": payload(detectors.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged full parent gives the actual canonical tree tensor wave operator on the fixed CD clock. A nonempty but restricted compact transverse-traceless detector class has two exact homogeneous-mode moment compatibilities, so its advanced stress test vanishes near both time endpoints and annihilates every free tensor initial-data contribution. A continuous all-momentum wave energy estimate, including zero momentum, bounds the adjoint stress-test norm by5000 times the declared detector norm. Combining it with the actual all-order CD stress covariance gives source-only tree-propagated canonical tensor standard deviation below5e-372 D(q), and dimensionless metric standard deviation below1e-771 D(q). No factorized zero quantum metric state, arbitrary gravity switch, full tensor self-energy control or finite-coupling remainder is asserted. Original V/G/B and P8 remain OPEN.",
        "not_established": [
            "An unrestricted finite-initial-time quantum metric state, all intrinsic or mixed correlations, or a Proca/CD boundary-divergence theorem",
            "The complete tensor self-energy inverse, its finite-coupling remainder or an interacting scalar/mixed/nonlinear stochastic background",
            "Physical cutoff/threshold/heavy or omitted-order matching, full parent measure/loops or an interacting vacuum amplitude",
            "Finite-gravity IR/Regge and original V/G/B or P8 closure",
        ],
        "verification_boundary": "Exact full-parent tensor action and canonical source, polynomial/orthogonal TT projectors, weighted Green/Wronskian identities and complete energy/derivative majorants. Independent actual-wave ODE, compact detector, incompatible-boundary, rotated-polarization, energy and normalization fixtures supplement the written continuous argument. The bound is a restricted source-only leading component, not a full quantum state. Native, direct science, ordinary and CLI use original SymPy; only full regression uses the audited exact GCD adapter. Analytic and distributional arguments are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The actual compatible tensor-response report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.187.ACTUAL_CD_BOUNDARY_COMPATIBLE_TREE_TENSOR_PROPAGATED_STRESS_NOISE_AND_ALL_MOMENTUM_ADJOINT_BOUND replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
