"""Read-only restricted global Euclidean control of the actual one-loop vacuum."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_quantum_map import verify as parent

from . import audit, calibration, dyson, kernel, map_domain

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-one-loop-Euclidean-control.json"
PARENT_SHA = "2560f29e703b4034773b85dc9fd7fa6d0a21d16f983daf985f39ee984963fee2"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_euclidean_control/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen actual one-loop quantum field-map report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_118_fully_rebuilt": PARENT_SHA,
        "same_full_parent_on_shell_self_energy_and_fixed_counterterm": True,
        "same_actual_MSbar_composite_coordinate": True,
        "proper_on_shell_amplitude_transfer_not_replaced_by_global_coordinate_norm": True,
        "no_old_rolling_state_or_finite_gravity_transfer": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = {k: bool(v) for k, v in audit.gates().items()}
    if not all(v is True for v in gates.values()):
        raise ValueError("An actual restricted Euclidean-control gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.119.POLYNOMIAL_VACUUM_ONE_LOOP_EUCLIDEAN_CONTROL",
        "date": "2026-09-09",
        "status": "COMPLETE_ACTUAL_GLOBAL_EUCLIDEAN_ONE_LOOP_RELATIVE_KERNEL_BOUND_SHARP_FIXED_NORMALIZATION_RESTRICTED_INSERTION_TAILS_AND_EXPLICIT_COMPOSITE_COORDINATE_NONUNIFORMITY; NOT_FULL_HIGHER_LOOP_AMPLITUDE_REFLECTION_POSITIVITY_TIMELIKE_CONTOUR_PHYSICAL_CUTOFF_ROW_NO_GO_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/kernel.md",
            "notes/dyson.md",
            "notes/map_domain.md",
            "notes/bounds.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": len(rows),
        "proof_checks": gates,
        "actual_global_Euclidean_kernel": serialize(payload(kernel.data())),
        "parent_light_and_stress_normalization": serialize(
            kernel.data()["same_light_and_stress_normalization"]
        ),
        "restricted_known_kernel_insertion_tails": serialize(
            {
                "general": payload(dyson.data()),
                "retained_order_zero": calibration.insertion_tail(0),
                "retained_order_one": calibration.insertion_tail(1),
                "retained_order_four": calibration.insertion_tail(4),
            }
        ),
        "actual_composite_coordinate_domain": serialize(payload(map_domain.data())),
        "actual_rational_global_and_pointwise_bounds": serialize(
            {
                "global": payload(calibration.data()),
                "origin": calibration.point(0),
                "unit_Euclidean_invariant": calibration.point(1),
                "heavy_scale": calibration.point(
                    calibration.data()["actual_heavy_mass_squared"]
                ),
                "large_composite_witness": calibration.point(10**400),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete fixed one-loop self-energy has a positive increasing relative Euclidean multiplier, uniformly below 3 times 10^-208, with sharp supremum equal to the independently checked residue/stress normalization. Positive form bounds and explicit geometric tails control repetitions of that kernel only. A logarithmically divergent local quartic bubble prevents promoting this comparison to a full renormalized b2 bound. The actual derivative composite coordinate stays small through the heavy window but has an explicit normalization-specific large-momentum nonuniformity witness. Neither this witness nor Euclidean positivity determines a physical cutoff, UV no-go or original P8 closure.",
        "not_established": [
            "Primitive higher-loop self-energies or vertices, their counterterms or a full higher-loop forward-coefficient bound",
            "Reflection positivity, timelike high-energy analyticity, heavy-resonance control or an all-energy UV completion",
            "A normalization-independent physical cutoff or instability from the derivative composite coordinate",
            "Finite-gravity IR/Regge contour and Delta_grav",
            "Old rolling-state transfer, absolute coupled renormalization or a common controlled bounce parent",
            "Exclusion or completion of an entire original P8 row, full V/G/B or original P8 closure",
        ],
        "verification_boundary": "Native exact full-kernel continuation, positive primitive and derivative identities, sharp limit, both stress-triangle normalization, known-kernel geometric coefficients/tails, ultraviolet-divergence control, actual composite completed square and rational bounds. Continuous dominated-limit, Fourier-multiplier and source-domain arguments are explicit source-pinned written proofs, not proof-assistant formalized, nonperturbative or peer reviewed. Every ancestor, own source and report field is rebuilt read-only with unmodified scientific SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The restricted global Euclidean one-loop report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.119.POLYNOMIAL_VACUUM_ONE_LOOP_EUCLIDEAN_CONTROL replay passed; full V/G/B and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
