"""Read-only regular analytic classical affine vacuum-domain certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import connection
from p8_affine import verify as certificate
from p8_exceptional_vacuum import verify as parent

from . import audit, bounds, calibration, family, local, lower

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "regular-analytic-affine-vacuum-lift.json"
PARENT_SHA = "2de0cb384f0c6361a106802a35934ed30572f015ce230fed07c7085298fc4355"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_affine_vacuum_domain/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {key: value for key, value in d.items() if key != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen exceptional vacuum tree-matching report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_108_fully_rebuilt": PARENT_SHA,
        "original_retuned_classical_action_and_full_literal_affine_dictionary_rebuilt": True,
        "new_analytic_family_and_auxiliary_boundary_do_not_edit_old_actions": True,
        "no_prior_kinetic_affine_parent_or_quantum_counterterm_transfer": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = {key: bool(value) for key, value in audit.gates().items()}
    if not all(value is True for value in gates.values()):
        raise ValueError("A regular affine vacuum-domain proof gate failed")
    d = family.data()
    R = d["R"]
    X = family.X
    return {
        "schema": 1,
        "claim": "P8-S6.109.REGULAR_ANALYTIC_CLASSICAL_AFFINE_VACUUM_LIFT",
        "date": "2026-09-09",
        "status": "COMPLETE_NEW_RANK_REGULAR_ANALYTIC_CLASSICAL_REDUCED_AND_AUXILIARY_AFFINE_VACUUM_TO_CLOCK_DOMAIN_FULL_PROJECTIVE_QUOTIENT_REGULAR_LOWER_INTEGRAL_FINITE_CLOCK_JETS_AND_SAME_MASSIVE_TREE_CONTACT; NOT_UNCHANGED_KINETIC_AFFINE_PARENT_PROPAGATING_UV_PARENT_QUANTUM_MEASURE_FULL_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/domain.md",
            "notes/norms.md",
            "notes/lift.md",
            "notes/lower.md",
            "notes/covariance.md",
            "notes/vacuum.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": sum(
            v.rows * v.cols if isinstance(v, sp.MatrixBase) else 1
            for v in rows.values()
        ),
        "proof_checks": gates,
        "new_analytic_reduced_family": serialize(d),
        "old_affine_extension_obstructions": serialize(
            payload(family.old_obstruction())
        ),
        "actual_principal_affine_lift": serialize(
            {
                "p": sp.sqrt(R) / 2,
                "c": -(sp.sqrt(R) - R) / X,
                "F4": (1 - R) / (2 * X * X),
                "Delta": sp.Integer(1),
                "additional_quotient_factor": 2 * R - 1,
                "physical_matter_metric_unchanged": True,
                "Proca_is_standalone_connection_independent": True,
            }
        ),
        "regular_lower_dictionary": serialize(payload(lower.data())),
        "general_gradient_quotient": serialize(
            {
                "unrestricted_components": 64,
                "projective_gauge_dimension": 4,
                "quotient_dimension": 60,
                "determinant": connection.quotient()["determinant"],
                "independent_full_matrix_cases": [
                    "timelike generic rational calibration",
                    "spacelike generic rational calibration",
                    "actual nonzero null u=0 n=1024 removable coefficients",
                    "actual zero-gradient vacuum coefficients",
                ],
                "continuous_domain_proof": "Lorentz-unimodular quotient covariance plus analytic identity on star-shaped gradient domain",
                "no_propagating_heavy_spectrum_inferred": True,
            }
        ),
        "uniform_clock_tube_and_domain_bounds": serialize(
            {
                "switch": payload(bounds.switch_jet_bounds()),
                "tube": bounds.tube(),
                "connected_domain": bounds.domain(),
            }
        ),
        "actual_vacuum_jets": serialize(payload(local.data())),
        "exact_calibrations": serialize(
            {
                "clock": calibration.point(),
                "vacuum": calibration.point(1024, 0),
                "spacelike": calibration.point(1024, -sp.Rational(1, 8192)),
                "next_even_order": calibration.point(1026, 0),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "A new rational-step/Gaussian analytic light-action family retains the actual clock jets and massive vacuum contact with uniform weighted four-jet tube control. Its full auxiliary connection quotient stays invertible throughout a connected domain, including null and zero gradients. A vacuum-regular coefficient integral supplies the complete lower dictionary and exact classical stationary-action equivalence. The direct older extension's rank-loss and clock-boundary obstructions are exhibited rather than hidden. The changed auxiliary parent is not the frozen kinetic completion, a propagating heavy UV parent or a full V/G/B result.",
        "not_established": [
            "Identification of the changed auxiliary boundary with the unchanged old kinetic affine completion or the separate healthy massive vacuum exchange model",
            "Propagating heavy masses, signs or a UV spectrum, interacting dispersion/cut control, finite-gravity Regge/IR-contour errors or complete V/G/B",
            "Uniform-in-n full auxiliary gap, every boosted coordinate inverse norm, off-clock scalar stability or nonlinear/quantum continuation",
            "Quantum path-integral measure equivalence, renormalized absolute sources, old state/counterterm reassignment, loop or omitted-operator bounds",
            "Original P8(b) classification or original P8 closure",
        ],
        "verification_boundary": "Native generic principal and complete lower affine identities, independent literal general-gradient 64-component Hessians and full quotient determinants, actual local scalar/affine jets and exact derivative/envelope recurrences. Continuous domain, Lorentz covariance, real-analytic extension, parameter integral and classical stationary-variation arguments are written/source-pinned, not proof-assistant formalized or peer reviewed. Every ancestor, own source/proof/test byte and report field is read-only rebuilt with unmodified scientific SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The regular analytic affine vacuum-lift report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.109.REGULAR_ANALYTIC_CLASSICAL_AFFINE_VACUUM_LIFT replay passed; propagating UV parent, full V/G/B and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
