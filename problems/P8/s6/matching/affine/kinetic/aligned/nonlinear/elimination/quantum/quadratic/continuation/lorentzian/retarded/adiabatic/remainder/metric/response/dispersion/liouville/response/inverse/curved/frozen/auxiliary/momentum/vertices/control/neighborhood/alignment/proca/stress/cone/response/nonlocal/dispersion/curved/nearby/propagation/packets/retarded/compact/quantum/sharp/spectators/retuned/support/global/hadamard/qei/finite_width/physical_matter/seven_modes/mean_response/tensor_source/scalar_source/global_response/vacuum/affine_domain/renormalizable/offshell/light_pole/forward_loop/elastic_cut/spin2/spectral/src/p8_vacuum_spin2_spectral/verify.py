"""Read-only actual spin-two spectral cuts and slope-moment bounds."""

import argparse
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_spin2 import verify as parent

from . import audit, calibration, cut, dispersion, moments, partial_wave

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-spin-two-spectral-cuts.json"
PARENT_SHA = "7bcd721b92f1711af646b4810c451f7c42d408fb2ec52b110f6792587f25aa77"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_spin2_spectral/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen actual spin-two form-factor report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_115_fully_rebuilt": PARENT_SHA,
        "same_actual_two_triangle_form_factor_and_once_fixed_Ward_counterterm": True,
        "same_light_mass_one_and_polynomial_cubic_coupling": True,
        "no_new_matching_action_or_Regge_ansatz": True,
        "heavy_decay_not_replaced_by_exact_stable_asymptotic_particle": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = {k: bool(v) for k, v in audit.gates().items()}
    if not all(v is True for v in gates.values()):
        raise ValueError("An actual spin-two spectral-cut gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.116.POLYNOMIAL_VACUUM_SPIN_TWO_SPECTRAL_CUTS",
        "date": "2026-09-09",
        "status": "COMPLETE_ACTUAL_ONE_LOOP_POSITIVE_VERTEX_SPECTRAL_REPRESENTATION_BOTH_UNEQUAL_MASS_CUTS_INDEPENDENT_SPIN_TWO_PHASE_NORMALIZATION_AND_STRICT_LOW_WINDOW_FULL_SLOPE_COMPARISON; NOT_EXACT_STABLE_HEAVY_SPECTRUM_FULL_GRAVITY_REGGE_DELTA_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/dispersion.md",
            "notes/cuts.md",
            "notes/moments.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": len(rows),
        "proof_checks": gates,
        "actual_positive_vertex_spectral_representation": serialize(
            payload(dispersion.data())
        ),
        "both_actual_unequal_mass_cut_densities": serialize(payload(cut.data())),
        "independent_spin_two_phase_space_and_positive_series": serialize(
            payload(partial_wave.data())
        ),
        "strict_low_window_and_full_species_moments": serialize(
            payload(moments.data())
        ),
        "exact_stable_spectral_calibrations": serialize(
            {
                "threshold": calibration.point(4),
                "interior": calibration.point(5),
                "upper_endpoint": calibration.point(6),
                "rational_interior": calibration.point(Fraction(9, 2)),
                "Q2_two": calibration.q2_enclosure(2, 3),
                "Q2_ten": calibration.q2_enclosure(10, 3),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The actual one-loop spin-two vertex has an exact positive once-subtracted spectral representation with distinct light and zeroth-order heavy stress cuts. Both densities follow from anchored triangle primitives and independently normalized integer-spin two-particle projections. Positive rational series enclosures avoid large-mass logarithmic cancellation. The cut from transfer four to six is strictly positive but below 10^-196 of the total slope; each full species moment exceeds one twentieth of the total. The heavy threshold is fixed-order, not an exact stable-particle assertion. No full scattering, Regge bound, Delta_grav or V/G/B closure is inferred.",
        "not_established": [
            "All-orders vertex spectral measure, exact stable-heavy threshold or uniform higher-loop threshold approximation",
            "Full gravitational scattering amplitude, massless cuts or graviton-loop error bounds",
            "Complex-spin continuation, Regge residue and trajectory, high-energy contour or numerical Delta_grav",
            "Quantum derivative-field-map equivalence or a common controlled rolling-bounce parent",
            "Exclusion or completion of an entire original P8 ladder row",
            "Full V/G/B, original P8(b) classification or original P8 closure",
        ],
        "verification_boundary": "Native exact pushforward identities, both actual mass assignments, anchored cut and angular primitives, phase-space factors, positive-series coefficients, parameter and transfer moments and actual rational hierarchy bounds. Continuous spectral, branch and perturbative-cut arguments are explicit source-pinned written proofs, not proof-assistant formalized or peer reviewed. Every ancestor, own source and report field is read-only rebuilt with unmodified scientific SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The polynomial-vacuum spin-two spectral report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.116.POLYNOMIAL_VACUUM_SPIN_TWO_SPECTRAL_CUTS replay passed; full V/G/B and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
