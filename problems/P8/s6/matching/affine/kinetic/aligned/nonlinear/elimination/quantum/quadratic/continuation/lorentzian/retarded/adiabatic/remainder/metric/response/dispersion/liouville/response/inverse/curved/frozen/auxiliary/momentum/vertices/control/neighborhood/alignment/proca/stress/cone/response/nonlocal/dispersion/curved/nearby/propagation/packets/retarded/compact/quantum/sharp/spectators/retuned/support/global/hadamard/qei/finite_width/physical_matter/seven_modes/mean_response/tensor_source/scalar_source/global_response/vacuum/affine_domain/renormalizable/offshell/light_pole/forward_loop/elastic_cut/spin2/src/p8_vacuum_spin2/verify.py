"""Read-only actual spin-two matter form factor and partial graviton-pole correction."""

import argparse
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_elastic_cut import verify as parent

from . import audit, calibration, pole, stress, triangles

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-spin-two-form-factor.json"
PARENT_SHA = "558e10ebc33b3a4ff9a0b09a9acf479cd18e1d6f2c2b62d014c0c0b855701f12"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_spin2/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {key: value for key, value in d.items() if key != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen complete one-loop elastic-cut report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_114_fully_rebuilt": PARENT_SHA,
        "actual_complete_one_loop_real_coefficient_and_elastic_cut_retained": True,
        "same_polynomial_model_light_mass_one_and_unit_residue": True,
        "same_once_fixed_light_kinetic_counterterm_for_gravitational_charge": True,
        "specified_minimal_Einstein_extension_not_affine_bounce_transfer": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = {key: bool(value) for key, value in audit.gates().items()}
    if not all(value is True for value in gates.values()):
        raise ValueError("An actual spin-two matter-form-factor gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.115.POLYNOMIAL_VACUUM_SPIN_TWO_FORM_FACTOR",
        "date": "2026-09-09",
        "status": "COMPLETE_ACTUAL_ONE_MATTER_LOOP_SPIN_TWO_FORM_FACTOR_BOTH_MASS_INSERTIONS_SAME_WARD_COUNTERTERM_POSITIVE_SLOPE_AND_NEGATIVE_FINITE_T_CHANNEL_VERTEX_CORRECTION_WITH_TRANSFER_REMAINDER; NOT_FULL_GRAVITY_AMPLITUDE_REGGE_DELTA_COMMON_BOUNCE_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/stress.md",
            "notes/triangles.md",
            "notes/pole.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": len(rows),
        "proof_checks": gates,
        "literal_canonical_stress_projection": serialize(payload(stress.data())),
        "both_actual_triangles_Ward_identity_slope_and_transfer_domain": serialize(
            payload(triangles.data())
        ),
        "literal_graviton_contraction_and_fixed_negative_t_limit": serialize(
            payload(pole.data())
        ),
        "exact_selected_and_transfer_calibrations": serialize(
            {
                "actual": calibration.point(),
                "smaller_transfer_disc": calibration.point(10**400, Fraction(1, 2)),
                "larger_controlled_transfer_disc": calibration.point(10**400, 3),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The actual polynomial vacuum with the specified minimal Einstein extension has a one-loop spin-two form factor from both mixed stress triangles. Its zero-transfer value is normalized by the already fixed light-residue counterterm; its positive slope is below 10^-405. After fixed-negative-transfer tree-pole subtraction, the t-channel matter vertex supplies a negative finite b2 correction of magnitude below 10^-1205 at Planck mass 10^400, with an explicit transfer-limit remainder. This is only that diagram class and order, not the full gravitational coefficient, a Regge bound, Delta_grav or V/G/B closure.",
        "not_established": [
            "Complete gravitational scattering amplitude or all exchange-channel coefficients at the retained inverse-Planck order",
            "Bounds on pure graviton loops, massless multiparticle cuts, higher matter loops or additional curvature matching",
            "High-energy Regge residue and trajectory, contour remainder or a numerical Delta_grav",
            "Quantum field-map equivalence or a common finite-gravity rolling-bounce parent",
            "Exclusion or completion of an entire original P8 ladder row",
            "Full V/G/B, original P8(b) classification or original P8 closure",
        ],
        "verification_boundary": "Native exact canonical stress, complete triangle routing and mass assignments, actual parent Ward and self-energy kernel checks, parameter integrals, tensor contraction, loop-order separation and rational transfer bounds. Covariant loop reduction, restricted diagram completeness and continuous inequalities are explicit primary-source-pinned written proofs, not proof-assistant formalized or peer reviewed. Every ancestor, own source and report field is read-only rebuilt with unmodified scientific SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The polynomial-vacuum spin-two report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.115.POLYNOMIAL_VACUUM_SPIN_TWO_FORM_FACTOR replay passed; full V/G/B and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
