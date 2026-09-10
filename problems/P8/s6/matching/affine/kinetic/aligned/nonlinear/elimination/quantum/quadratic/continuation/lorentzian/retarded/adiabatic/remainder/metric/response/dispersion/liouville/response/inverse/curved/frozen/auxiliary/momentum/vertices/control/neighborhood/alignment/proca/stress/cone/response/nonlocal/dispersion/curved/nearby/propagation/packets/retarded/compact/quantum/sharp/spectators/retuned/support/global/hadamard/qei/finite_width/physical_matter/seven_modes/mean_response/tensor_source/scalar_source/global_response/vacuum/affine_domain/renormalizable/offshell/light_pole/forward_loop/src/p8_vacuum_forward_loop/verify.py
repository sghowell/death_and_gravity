"""Read-only complete one-loop polynomial-vacuum forward-coefficient certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as certificate
from p8_vacuum_light_pole import verify as parent

from . import angular, audit, calibration, domain, leading, normalization, subtraction

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-one-loop-forward-coefficient.json"
PARENT_SHA = "a4c09a394e9fa8a00aa3338c786ff68f9ba9780c1868fbb2d21bbb8cbeaa224d"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_forward_loop/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {key: value for key, value in d.items() if key != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen actual one-loop light-pole report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_112_fully_rebuilt": PARENT_SHA,
        "same_S6_110_full_polynomial_model_and_Hessian_rebuilt": True,
        "actual_light_on_shell_mass_one_and_unit_residue_scheme_retained": True,
        "classical_field_map_not_used_to_transfer_quantum_loops_or_measure": True,
        "no_old_cosmological_state_or_counterterm_transfer": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = {key: bool(value) for key, value in audit.gates().items()}
    if not all(value is True for value in gates.values()):
        raise ValueError("An actual complete one-loop forward-coefficient gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.113.POLYNOMIAL_VACUUM_ONE_LOOP_FORWARD_COEFFICIENT",
        "date": "2026-09-09",
        "status": "COMPLETE_ACTUAL_SAME_MODEL_FULL_ONE_LOOP_FORWARD_B2_ERROR_BELOW_ONE_MILLIONTH_OF_POSITIVE_EXACT_TREE_VALUE_WITH_FULL_VERTICES_ANALYTIC_DISC_UNBOUNDED_RADIAL_INTEGRALS_AND_FIXED_LOCAL_COUNTERTERMS; NOT_HIGHER_LOOP_ALL_ENERGY_CONTOUR_FINITE_GRAVITY_COMMON_BOUNCE_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/normalization.md",
            "notes/analytic.md",
            "notes/angular.md",
            "notes/radial.md",
            "notes/subtraction.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": sum(
            value.rows * value.cols if isinstance(value, sp.MatrixBase) else 1
            for value in rows.values()
        ),
        "proof_checks": gates,
        "actual_full_vertex_normalization": serialize(payload(normalization.data())),
        "holomorphic_forward_domain": serialize(payload(domain.data())),
        "uniform_angular_remainder": serialize(payload(angular.data())),
        "subtracted_radial_derivative_bound": serialize(payload(leading.data())),
        "fixed_local_counterterm_scheme": serialize(payload(subtraction.data())),
        "exact_error_calibrations": serialize(
            {
                "actual_parameters_and_inequalities": payload(calibration.data()),
                "one_millionth": calibration.point(),
                "looser": calibration.point(sp.Rational(1, 10**5)),
                "unproved_stricter_tolerance": calibration.point(sp.Rational(1, 10**8)),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete one-loop forward coefficient of the specified polynomial vacuum is bounded in one fixed local full-model subtraction scheme. All three channels and exact nonlocal heavy vertices remain inside the unrestricted loop integral. A full-parameter analytic disc, uniform angular Cauchy bound and explicitly UV-subtracted radial derivative estimate give an exact error below 10^-6 times the actual positive tree b2=4lambda. The tree-plus-one-loop coefficient is therefore positive. No all-orders error, high-energy contour, finite-gravity or common-bounce matching is established.",
        "not_established": [
            "All higher-loop or nonperturbative errors, convergence or an all-energy UV completion",
            "Complete controlled absorptive density through the unstable heavy resonance or complex high-energy dispersion contour",
            "Finite-gravity infrared/Regge remainder or the G test",
            "Quantum field-map equivalence, corrected cosmological cones, absolute coupled sources or the common-parent B test",
            "Exclusion or completion of an entire original P8 ladder row",
            "Full V/G/B, original P8(b) classification or original P8 closure",
        ],
        "verification_boundary": "Native full-Hessian quartic tensor and trace differentiation, exact channel counterterm cancellation, forward kinematics, untruncated external-shift identities, convergent radial integrals and actual rational error margins. Continuum analyticity, dominated integration, angular averaging and Cauchy arguments are explicit source-pinned written proofs, not proof-assistant formalized or peer reviewed. Every ancestor, own source and report field is read-only rebuilt with unmodified scientific SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The polynomial-vacuum forward-loop report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.113.POLYNOMIAL_VACUUM_ONE_LOOP_FORWARD_COEFFICIENT replay passed; full V/G/B and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
