"""Read-only global finite relative-Proca first-order mean-response certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as certificate
from p8_proca_seven_modes import verify as parent

from . import anchor, audit, band, mean, stress, tails

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "global-relative-Proca-mean-response.json"
PARENT_SHA = "498312d2df8964c169851c242eb9fb898edd3ea763d0efb52846fac02a9f69ca"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_proca_mean_response/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen actual seven-mode state certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_103_fully_rebuilt": PARENT_SHA,
        "same_original_global_retuned_action_and_physical_metric": True,
        "same_scalar_tensor_factors_and_counterterms_in_both_states": True,
        "finite_Proca_state_difference_not_a_reassigned_reference_tadpole": True,
    }


def payload(data):
    return {k: v for k, v in data.items() if k != "checks"}


@cache
def build_report():
    prior = prior_checks()
    identities = audit.residuals()
    gates = {key: bool(value) for key, value in audit.gates().items()}
    if not all(value is True for value in gates.values()):
        raise ValueError("A global relative-Proca mean-response proof gate failed")
    exact = certificate.certify_residuals(identities)
    return {
        "schema": 1,
        "claim": "P8-S6.104.GLOBAL_RELATIVE_PROCA_MEAN_RESPONSE",
        "date": "2026-09-09",
        "status": "COMPLETE_FINITE_HADAMARD_PROCA_STATE_DIFFERENCE_CONSERVED_PHYSICAL_STRESS_AND_ACTUAL_GLOBAL_BOUNDED_LEADING_HOMOGENEOUS_MEAN_RESPONSE_WITH_MATTER_FIELD_CONNECTION_VARIATION_AND_EXACT_FRAME_COMPLETE_BOUNCE_REPRESENTATIVE; NOT_AN_ABSOLUTE_OR_EXACT_SEMICLASSICAL_SOLUTION_HIGHER_ORDER_CONTROL_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/state.md",
            "notes/stress.md",
            "notes/mean.md",
            "notes/calibration.md",
            "notes/tails.md",
            "notes/representative.md",
            "notes/scope.md",
        ],
        "exact_residuals": exact,
        "named_exact_check_count": len(identities),
        "checked_scalar_entries": sum(
            v.rows * v.cols if isinstance(v, sp.MatrixBase) else 1
            for v in identities.values()
        ),
        "proof_checks": gates,
        "physical_relative_Proca_stress": serialize(payload(stress.data())),
        "actual_leading_mean_response": serialize(
            {**payload(mean.data()), "center_control": payload(mean.center_control())}
        ),
        "physical_Ward_identity": serialize(payload(mean.ward())),
        "fixed_band_and_interval_bounds": serialize(
            {**payload(band.data()), "explicit_calibration": band.calibrated()}
        ),
        "global_weighted_response_and_frame": serialize(payload(tails.data())),
        "independent_center_stress_jets": serialize(payload(anchor.data())),
        "controls": serialize(audit.controls()),
        "verdict": "A fixed positive compact-band addition to the actual Proca covariance gives a smooth Hadamard state difference, finite conserved physical stress and a uniquely determined leading homogeneous metric/matter response of the original retuned action. The actual spatial-only lapse chain retains pressure, the full six-channel induced matter Hessian is independently recovered, and the matter-field response includes the connection variation required by the homogeneous Ward identity. Independent compact-band and exact center-jet bounds are supplemented by global weighted estimates and both time tails. The exact-frame first-order profile has a complete bouncing metric representative, but neither that metric nor a newly co-evolved quantum state is claimed to solve an absolute semiclassical equation. The common vacuum source, nonlinear/higher-loop errors and V/G/B matching remain open.",
        "not_established": [
            "The common absolute vacuum tadpole, full seven-mode covariant stress renormalization or interacting quantum Ward identities",
            "An actual finite-amplitude state-difference or absolute semiclassical solution, a quantum state co-evolved on the representative, or a nonlinear/higher-loop residual bound",
            "Finite total energy of the homogeneous state on infinite R3, an arbitrary-amplitude perturbative regime, or a physical-matter QEI",
            "A heavy-spectrum or finite-frequency EFT cutoff, omitted-operator/loop error control, actual vacuum/finite-gravity/common-parent V/G/B matching, or UV completion",
            "Original P8(b) classification or original P8 closure",
        ],
        "verification_boundary": "Native exact constrained Proca stress and conservation, correct spatial-only lapse variation with a nonzero four-dimensional-rescaling control, full six-channel physical-density bridge, literal retuned-action mean/charge coefficients, complete homogeneous relative Ward identity, independent center stress jets and proper-time lapse factors, exact rational whole-interval/global integral bounds and negative-input controls. Written state, differential-equation, tail and geometric representative proofs are source-pinned but not proof-assistant formalized. They certify leading state-difference equations, not uncomputed common or higher-order sources. Every ancestor, own source/proof/test byte and report field is independently replayed using native scientific SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The global relative-Proca mean-response report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.104.GLOBAL_RELATIVE_PROCA_MEAN_RESPONSE replay passed; absolute source, V/G/B and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
