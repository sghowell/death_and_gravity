"""Read-only actual off-shell quartic vacuum tree-matching certificate."""

import argparse
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_polynomial_vacuum import verify as parent

from . import audit, band, calibration, jets, matching, norms

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "off-shell-vacuum-quartic-matching.json"
PARENT_SHA = "688e57e0907d1311c9e30c06966ad188d29ad15a0e72f4452b3e1c4440ef0230"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_offshell_vacuum/*.py"))
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
            "The frozen polynomial vacuum quantum-preparation report changed"
        )
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_110_fully_rebuilt": PARENT_SHA,
        "actual_S6_109_canonical_affine_vacuum_quartic_rebuilt": True,
        "S6_108_fixed_interaction_limit_and_full_tree_exchange_rebuilt": True,
        "no_old_state_counterterm_or_finite_M_parent_transfer": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = {key: bool(value) for key, value in audit.gates().items()}
    if not all(value is True for value in gates.values()):
        raise ValueError("An actual off-shell vacuum tree-matching gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.111.OFF_SHELL_VACUUM_QUARTIC_MATCHING",
        "date": "2026-09-09",
        "status": "COMPLETE_ACTUAL_OFF_SHELL_QUARTIC_VACUUM_TREE_MATCH_EXPLICIT_CUBIC_FIELD_MAP_FULL_BOUNDARY_CURRENT_4D_JET_IDENTITY_AND_CONTROLLED_COMMON_SCHWARTZ_ACTION_REMAINDER; NOT_ALL_VERTEX_QUANTUM_FINITE_GRAVITY_OR_COMMON_BOUNCE_PARENT_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/identity.md",
            "notes/field_map.md",
            "notes/norms.md",
            "notes/common_class.md",
            "notes/scope.md",
            "notes/literature.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": len(rows),
        "proof_checks": gates,
        "actual_affine_vacuum_target": serialize(payload(matching.actual_target())),
        "off_shell_pointwise_identity": serialize(payload(jets.data())),
        "full_centered_resolvent": serialize(payload(matching.data())),
        "local_jet_and_substitution_norms": serialize(payload(norms.data())),
        "common_Schwartz_domain": serialize(payload(band.data())),
        "exact_calibrations": serialize(
            {
                "wide_tree_window": calibration.point(),
                "common_class_window": calibration.point(38),
                "half_jet_cube": calibration.point(38, Fraction(1, 2)),
                "zero_jet_bound": calibration.point(1, 0),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The actual canonical affine vacuum quartic matches the polynomial heavy model off shell through an explicit cubic local derivative field map and full divergence. The complete four-dimensional sixth-jet identity is exact. Higher field powers and the untruncated heavy inverse have separate quantitative bounds and a nonempty common Schwartz/Fourier domain yields an integrated tree-action error below 10^-800 times the squared light norm. No global phase-space map, all-vertex or quantum matching, finite-gravity contour or full rolling-bounce parent is identified.",
        "not_established": [
            "Equality of all higher vacuum vertices, the entire finite-M analytic affine action or the earlier exponential kinetic parent",
            "Global nonlinear solution/phase-space invertibility of the derivative-dependent map or its use on the entire rolling clock",
            "Automatic Jacobian, source, state or counterterm transfer, complete renormalized scattering loops or quantum pole/residue data",
            "Holomorphic loop extension of a tree spectral window or a restriction of actual loop momenta to the classical test class",
            "Finite-gravity IR/Regge remainder, corrected cones, absolute cosmological sources or a controlled common propagating bounce parent",
            "Original P8(b) classification or original P8 closure",
        ],
        "verification_boundary": "Native full 4D independent field-jet differentiation and pointwise polynomial/current identities, actual affine vacuum coefficients, exact centered resolvent and rational coefficient majorants. Finite-polynomial norm induction, Fourier support, Schwartz integration and formal field-map arguments are written/source-pinned, not proof-assistant formalized or peer reviewed. Every ancestor, own source and report field is read-only rebuilt with unmodified scientific SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The off-shell vacuum quartic-matching report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.111.OFF_SHELL_VACUUM_QUARTIC_MATCHING replay passed; full V/G/B and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
