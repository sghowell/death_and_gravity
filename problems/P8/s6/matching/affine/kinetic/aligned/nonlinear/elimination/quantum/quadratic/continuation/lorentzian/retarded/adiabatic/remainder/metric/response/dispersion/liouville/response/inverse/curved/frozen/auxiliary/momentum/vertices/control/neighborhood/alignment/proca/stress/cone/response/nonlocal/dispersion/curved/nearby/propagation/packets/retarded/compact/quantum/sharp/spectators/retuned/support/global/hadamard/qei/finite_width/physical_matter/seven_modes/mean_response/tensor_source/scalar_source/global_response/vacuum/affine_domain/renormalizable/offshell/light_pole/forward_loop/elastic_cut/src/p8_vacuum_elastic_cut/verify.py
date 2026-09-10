"""Read-only actual one-loop elastic cut and finite-window improved coefficient."""

import argparse
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_forward_loop import verify as parent

from . import amplitude, angular, audit, calibration, cut, normalization

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-one-loop-elastic-cut.json"
PARENT_SHA = "336e83ee3053018d9ce8a84a90651e0ceab79eb2d176003c5e466faf126fc2c2"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_elastic_cut/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {key: value for key, value in d.items() if key != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError(
            "The frozen complete one-loop forward-coefficient report changed"
        )
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_113_fully_rebuilt": PARENT_SHA,
        "actual_complete_one_loop_real_coefficient_error_retained": True,
        "same_polynomial_model_light_mass_one_and_unit_residue": True,
        "same_fixed_real_local_counterterm_scheme": True,
        "no_old_cosmological_state_or_counterterm_transfer": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = {key: bool(value) for key, value in audit.gates().items()}
    if not all(value is True for value in gates.values()):
        raise ValueError("An actual one-loop elastic-cut gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.114.POLYNOMIAL_VACUUM_ONE_LOOP_ELASTIC_CUT",
        "date": "2026-09-09",
        "status": "COMPLETE_ACTUAL_ONE_LOOP_ELASTIC_LOW_ENERGY_CUT_EXACT_FULL_VERTEX_ANGULAR_INTEGRAL_IDENTICAL_PARTICLE_NORMALIZATION_STRICT_POSITIVE_WEIGHT_AND_POSITIVE_CUT_SUBTRACTED_TREE_PLUS_ONE_LOOP_COEFFICIENT; NOT_HIGHER_LOOP_GLOBAL_DISPERSION_FINITE_GRAVITY_COMMON_BOUNCE_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/amplitude.md",
            "notes/normalization.md",
            "notes/cut.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": len(rows),
        "proof_checks": gates,
        "actual_physical_tree_amplitude": serialize(payload(amplitude.data())),
        "exact_angular_integral": serialize(payload(angular.data())),
        "independent_cut_normalization": serialize(payload(normalization.data())),
        "actual_low_energy_cut_and_subtracted_margin": serialize(payload(cut.data())),
        "exact_invariant_calibrations": serialize(
            {
                "threshold": calibration.point(4),
                "interior": calibration.point(5),
                "upper_endpoint": calibration.point(6),
                "rational_interior": calibration.point(Fraction(9, 2)),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The actual full nonlocal tree vertex gives a strictly positive one-loop elastic discontinuity above the mass-one two-particle threshold and below invariant six. A closed anchored angular integral and independent identical-particle/bubble normalization retain every interference term. The computed low-energy forward cut lies between lambda squared over twenty and three lambda squared. Subtracting it together with the complete one-loop real-coefficient error leaves a strictly positive coefficient with the same one-millionth relative margin. This is finite-order evidence, not a global improved positivity or V/G/B closure verdict.",
        "not_established": [
            "All higher-loop corrections to the cut, real coefficient or scattering amplitude",
            "Global twice-subtracted dispersion equality, high-energy contour control or an all-energy UV completion",
            "A controlled spectral density through an unstable heavy resonance",
            "Finite-gravity infrared/Regge remainder, quantum field-map equivalence or a common rolling-bounce parent",
            "Exclusion or completion of an entire original P8 ladder row",
            "Full V/G/B, original P8(b) classification or original P8 closure",
        ],
        "verification_boundary": "Native exact physical tree amplitude, monotonicity identities, anchored angular primitives, phase-space and bubble-cut factors, threshold limits, positive weight integrals and actual rational error margins. Continuous inequalities and perturbative unitarity/dispersion-scope arguments are explicit primary-source-pinned written proofs, not proof-assistant formalized or peer reviewed. Every ancestor, own source and report field is read-only rebuilt with unmodified scientific SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The polynomial-vacuum elastic-cut report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.114.POLYNOMIAL_VACUUM_ONE_LOOP_ELASTIC_CUT replay passed; full V/G/B and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
