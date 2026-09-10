"""Read-only actual polynomial-vacuum one-loop subtraction-reference flow."""

import argparse
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_spin2_spectral import verify as parent

from . import audit, calibration, flow, subtraction, tensors

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-one-loop-reference-flow.json"
PARENT_SHA = "12944ab16a9efa8c2841ec456a8aee5d4ea022e65a333de4f7120901f4b42dd9"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_polynomial_vacuum_running/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen actual spin-two spectral report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_116_fully_rebuilt": PARENT_SHA,
        "same_polynomial_vacuum_and_fixed_one_loop_full_model_counterterms": True,
        "same_initial_finite_potential_contact_and_light_on_shell_conditions": True,
        "no_new_Regge_ansatz_or_common_bounce_transfer": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = {k: bool(v) for k, v in audit.gates().items()}
    if not all(v is True for v in gates.values()):
        raise ValueError("An actual one-loop reference-flow gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.117.POLYNOMIAL_VACUUM_ONE_LOOP_REFERENCE_FLOW",
        "date": "2026-09-09",
        "status": "COMPLETE_ACTUAL_ONE_LOOP_FULL_TENSOR_BETA_NORMALIZATION_ANCHORED_SUBTRACTION_REFERENCE_COMPENSATION_CLOSED_POSITIVE_MARGIN_FLOW_AND_UNIFORM_DECLARED_HIERARCHY_ENCLOSURES; NOT_HIGHER_LOOP_QUANTUM_FLOW_PLANCKIAN_SCATTERING_REGGE_DELTA_COMMON_BOUNCE_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/tensors.md",
            "notes/subtraction.md",
            "notes/flow.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": len(rows),
        "proof_checks": gates,
        "literal_two_field_beta_tensors": serialize(payload(tensors.data())),
        "anchored_reference_change_and_full_channel_compensation": serialize(
            payload(subtraction.data())
        ),
        "closed_one_loop_flow_and_positive_margin": serialize(payload(flow.data())),
        "actual_uniform_and_pointwise_reference_enclosures": serialize(
            {
                "whole_window": payload(calibration.data()),
                "initial_reference": calibration.point(0),
                "interior_reference": calibration.point(1),
                "rational_reference": calibration.point(Fraction(9, 2)),
                "upper_enclosure_reference": calibration.point(9),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The actual full-model one-loop counterterms agree with a literal complete two-field scalar tensor calculation. An anchored change of reference bubble is compensated in every scattering channel through one loop, without retuning the finite potential contact or physical light mass. The closed one-loop beta flow preserves the positive completed-square quartic margin. Rational enclosures include all subtraction references from the light mass to the selected Planck scale; the relative margin increase is below 10^-5 and dimensionless coupling changes below 10^-202. This is a truncated reference-flow statement, not higher-loop control, a Planckian amplitude, a quantum Landau-pole theorem or V/G/B closure.",
        "not_established": [
            "Higher-loop beta coefficients, finite-scheme conversion remainder or exact quantum reference flow",
            "All-energy scattering, exact heavy spectral behavior or a nonperturbative UV completion",
            "Graviton-loop, massless-cut, higher-curvature or Regge-contour bounds and numerical Delta_grav",
            "Quantum derivative-field-map equivalence or a common controlled rolling-bounce parent",
            "Exclusion or completion of an entire original P8 ladder row",
            "Full V/G/B, original P8(b) classification or original P8 closure",
        ],
        "verification_boundary": "Native exact full scalar tensor contractions, existing-counterterm normalization, convergent reference-bubble primitive and anchor, all-channel compensation, closed beta ODE and composition, positive-margin identities and actual rational window bounds. Continuous flow and perturbative scheme arguments are explicit source-pinned written proofs, not proof-assistant formalized or peer reviewed. Every ancestor, own source and report field is read-only rebuilt with unmodified scientific SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The polynomial-vacuum one-loop reference-flow report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.117.POLYNOMIAL_VACUUM_ONE_LOOP_REFERENCE_FLOW replay passed; full V/G/B and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
