"""Read-only actual tensor clock-source and global leading mean-response certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as certificate
from p8_proca_mean_response import verify as parent

from . import audit, bounds, calibration, chart, operator, response, state

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "actual-tensor-clock-mean-response.json"
PARENT_SHA = "b27aec2f3eed846ac1201bbb6d5470b6d7d3c98d439e3da48f708af66466e68c"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_tensor_mean_source/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError(
            "The frozen global relative-Proca response certificate changed"
        )
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_104_fully_rebuilt": PARENT_SHA,
        "same_original_global_retuned_action_and_physical_metric": True,
        "same_reference_scalar_and_vector_factors_and_counterterms": True,
        "new_named_smooth_tensor_covariance_addition_only": True,
    }


def payload(data):
    return {k: v for k, v in data.items() if k != "checks"}


@cache
def build_report():
    prior = prior_checks()
    identities = audit.residuals()
    gates = {key: bool(value) for key, value in audit.gates().items()}
    if not all(value is True for value in gates.values()):
        raise ValueError("An actual tensor-source proof gate failed")
    exact = certificate.certify_residuals(identities)
    return {
        "schema": 1,
        "claim": "P8-S6.105.ACTUAL_TENSOR_CLOCK_MEAN_RESPONSE",
        "date": "2026-09-09",
        "status": "COMPLETE_ACTUAL_HOMOGENEOUS_TENSOR_STATE_DIFFERENCE_LAPSE_PRESSURE_AND_CLOCK_SOURCE_WITH_NONLINEAR_CENTER_CHART_BRIDGE_AND_GLOBAL_BOUNDED_LEADING_METRIC_MATTER_RESPONSE; NOT_FULL_COVARIANT_GRAVITON_STRESS_AN_ABSOLUTE_OR_EXACT_SEMICLASSICAL_SOLUTION_HIGHER_ORDER_CONTROL_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/action.md",
            "notes/state.md",
            "notes/chart.md",
            "notes/mean.md",
            "notes/clock.md",
            "notes/bounds.md",
            "notes/scope.md",
        ],
        "exact_residuals": exact,
        "named_exact_check_count": len(identities),
        "checked_scalar_entries": sum(
            v.rows * v.cols if isinstance(v, sp.MatrixBase) else 1
            for v in identities.values()
        ),
        "proof_checks": gates,
        "actual_tensor_action_and_lapse_source": serialize(
            {
                **payload(operator.data()),
                "homogeneous_clock_identity": payload(operator.clock()),
            }
        ),
        "positive_actual_tensor_state_difference": serialize(
            {
                **payload(state.data()),
                "direct_clock_variation": payload(state.variational_clock()),
            }
        ),
        "actual_leading_tensor_mean_response": serialize(
            {
                **payload(response.data()),
                "independent_center_jets": payload(response.center()),
            }
        ),
        "actual_clock_and_joint_source_Ward_identity": serialize(
            payload(response.ward())
        ),
        "nonlinear_center_physical_density_chart": serialize(
            {**payload(chart.data()), "literal_metric_map": payload(chart.metric())}
        ),
        "global_response_and_clock_transfer_bounds": serialize(
            {**payload(bounds.data()), "explicit_calibration": calibration.calibrated()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "A fixed positive compact-band addition to both actual tensor covariances gives a smooth Hadamard state difference. Varying the original tensor coefficient before clock restriction yields the correct physical lapse, spatial pressure and nonzero clock source, not the minimally coupled Proca source. The complete nonlinear spatial chart reconciles the center matter-density observable with the earlier additive calculation. The same classical mean operator with new forcing has global weighted bounds, absolutely integrable clock exchange and a complete bouncing exact-frame representative. This is a leading homogeneous state-difference result; the common absolute source, covariant graviton renormalization, nonlinear and higher-loop errors and V/G/B matching remain open.",
        "not_established": [
            "A full covariant graviton or ghost/BRST stress renormalization, interacting Ward identity, coupled-scalar absolute source, or vanishing common vacuum tadpole",
            "An exact finite-amplitude or absolute semiclassical solution, a quantum state co-evolved on the representative, or nonlinear/higher-loop residual control",
            "An all-time additive canonical momentum chart bridge, finite total energy on infinite R3, a physical-matter QEI, or a UV no-go from a negative quadratic source",
            "A heavy-spectrum or finite-frequency cutoff, omitted-operator or loop bounds, actual vacuum/finite-gravity/common-parent V/G/B matching, or UV completion",
            "Original P8(b) classification or original P8 closure",
        ],
        "verification_boundary": "Native exact variation of the actual tensor coefficient and spatial-only frame, both TT polarizations, independent source/Hessian transport, induced matter connection and clock Ward identities, explicit nonlinear center metric/observable bridge, rational global and whole-line integral bounds, parity and negative-input controls. Written smooth-state, differential-equation and geometric arguments are source-pinned, not proof-assistant formalized. Every ancestor, own source/proof/test byte and report field is independently replayed with unmodified scientific SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The actual tensor clock-source report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.105.ACTUAL_TENSOR_CLOCK_MEAN_RESPONSE replay passed; absolute source, V/G/B and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
