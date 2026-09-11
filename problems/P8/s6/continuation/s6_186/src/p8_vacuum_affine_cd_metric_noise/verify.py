"""Read-only full local CD metric stress covariance and smearing bound."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_cd_source_noise import verify as parent
from p8_vacuum_affine_proca_gaussian import verify as gaussian

from . import audit, noise, oscillatory, reference, stress

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-cd-metric-noise.json"
PARENT_SHA = "39e536b6de14a434e8db2a69e65529723abbf8c5769629f5e1f3dc4019d97eb6"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_cd_metric_noise/*.py"))
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
        "S6_185_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("An actual CD metric-noise gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.186.ACTUAL_ALL_ORDER_CD_STATE_FULL_LOCAL_SPACETIME_SMEARED_PROCA_METRIC_STRESS_COVARIANCE_AND_NORMALIZATION_BOUND",
        "date": "2026-09-11",
        "status": "SAME_ALL_ORDER_CD_STATE_FULL_LOCAL_SMEARED_METRIC_STRESS_NOISE_BOUND; NOT_METRIC_SOLUTION_REDUCED_MIXED_MODE_STABILITY_QUANTUM_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/stress.md",
            "notes/reference.md",
            "notes/state.md",
            "notes/oscillatory.md",
            "notes/noise.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "full_stress_and_actual_state": serialize(
            {"stress": payload(stress.data()), "reference": payload(reference.data())}
        ),
        "complete_pair_smearing_and_noise": serialize(
            {"pair": payload(oscillatory.data()), "noise": payload(noise.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The actual unchanged all-order CD Gaussian vector state's full local metric stress has a continuous compact-spacetime smearing bound, retaining every tensor component, polarization, temporal constraint and internal/external momentum. A full analytic reference-phase estimate and the exact evolved-state remainder control the complete two-particle Wick kernel. Its standard deviation is below1e25 in the declared three-time-derivative and one-spatial-derivative tensor norm, below1e-775 for stress/kappa and below1e-375 for the Einstein-normalized metric test coupling. The same fixed mean prescription and scalar retuning are retained; c-numbers cancel only in centering. This is a conditional stress-noise input, not a fully reduced mixed-mode norm, metric solution, stability, physical cutoff, interacting quantum background or original V/G/B closure.",
        "not_established": [
            "A coincident stress variance, Gaussian higher stress statistics or a physical momentum cutoff",
            "Fully reduced canonical mixed-mode noise, a small or stable full inverse, spatial/nonlinear/stochastic metric solutions or finite-coupling remainders",
            "Complete interacting parent state/measure/loops and heavy/cutoff/threshold or omitted-order matching",
            "Complete vacuum cuts/contour/truncation, finite-gravity IR/Regge or original V/G/B closure",
        ],
        "verification_boundary": "Exact physical stress, constrained mode readouts, full integration-by-parts/Leibniz and radial identities; continuous complex coefficient and exact-state remainder bounds; complete two-particle covariance and smearing argument. Independent covariant metric variation, complex-time reference, actual preparation, polarization rotation, Fock-pair, radial and normalization fixtures supplement the proofs. Native, direct science, ordinary and CLI use original SymPy; only full regression uses the audited exact GCD adapter. Analytic and distributional arguments are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The actual CD metric-noise report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.186.ACTUAL_ALL_ORDER_CD_STATE_FULL_LOCAL_SPACETIME_SMEARED_PROCA_METRIC_STRESS_COVARIANCE_AND_NORMALIZATION_BOUND replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
