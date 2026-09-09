"""Read-only global scalar and seven-mode leading relative mean certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as certificate
from p8_scalar_mean_source import verify as parent

from . import audit, bounds, calibration, identities, joint, phase, source

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "global-seven-mode-relative-mean-control.json"
PARENT_SHA = "6a8deadd89d451245f8c3f1d3c9de93f6b85dcc1b165e40e8df826d0e71d6aa5"
sha, serialize = certificate.sha, certificate.serialize
rational_signature, matrix_signature = (
    parent.rational_signature,
    parent.matrix_signature,
)


def source_files():
    return (
        sorted(ROOT.glob("src/p8_scalar_global_response/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {key: value for key, value in d.items() if key != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen actual coupled scalar source certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_106_fully_rebuilt": PARENT_SHA,
        "same_actual_action_reference_and_common_counterterms": True,
        "same_scalar_covariance_family_on_a_smaller_global_amplitude_domain": True,
        "tensor_and_Proca_additions_and_sources_rebuilt_through_ancestry": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = {key: bool(value) for key, value in audit.gates().items()}
    if not all(value is True for value in gates.values()):
        raise ValueError("A global scalar or joint seven-mode proof gate failed")
    exact = certificate.certify_residuals(rows)
    ph = phase.data()
    src = source.data()
    return {
        "schema": 1,
        "claim": "P8-S6.107.GLOBAL_SEVEN_MODE_RELATIVE_MEAN_CONTROL",
        "date": "2026-09-09",
        "status": "COMPLETE_GLOBAL_ACTUAL_COUPLED_SCALAR_COVARIANCE_COMPENSATED_SAME_MEAN_RESPONSE_PHYSICAL_OBSERVABLE_AND_JOINT_SEVEN_MODE_FINITE_RELATIVE_LEADING_CONTROL_WITH_COMPLETE_BOUNCING_FRAME_REPRESENTATIVE; NOT_COMMON_ABSOLUTE_SOURCE_EXACT_SEE_HIGHER_ORDER_CORRECTED_CONES_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/phase.md",
            "notes/compensation.md",
            "notes/bounds.md",
            "notes/observable.md",
            "notes/joint.md",
            "notes/geometry.md",
            "notes/scope.md",
        ],
        "exact_residuals": exact,
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": sum(
            v.rows * v.cols if isinstance(v, sp.MatrixBase) else 1
            for v in rows.values()
        ),
        "proof_checks": gates,
        "actual_weighted_scalar_phase": serialize(
            {
                **{
                    key: value
                    for key, value in ph.items()
                    if not isinstance(value, sp.MatrixBase)
                },
                "complete_rational_matrix_fingerprints": {
                    key: matrix_signature(value)
                    for key, value in ph.items()
                    if isinstance(value, sp.MatrixBase)
                },
                "phase_weight": "W=diag(1,1,(2/5)t^-3,(2/5)t^-3) Z_old; t=1+u^2",
                "actual_CCR": "Omega_W=(2/5)Omega/t^3, not reset to canonical identity",
                "improved_weight": "(v,chi,sqrt(t)W_pv,sqrt(t)W_pchi)",
                "added_covariance": "same natural eta_S*b*I4 anchor and actual two-leg evolution",
            }
        ),
        "compensated_actual_mean_kernels": serialize(
            {
                "compensating_trace_coefficient": src["compensating_trace_coefficient"],
                "same_physical_trace": "d=delta_p+9H*<v^2>; delta_p=d-9H*<v^2>",
                "curvature_variance_Hessian": src["curvature_variance_Hessian"],
                "curvature_variance_derivative_Hessian": matrix_signature(
                    src["curvature_variance_derivative_Hessian"]
                ),
                "compensated_old_phase_sources": {
                    key: matrix_signature(v)
                    for key, v in src["compensated_old_phase_sources"].items()
                },
                "actual_weighted_response_kernels": {
                    key: matrix_signature(v)
                    for key, v in src["actual_weighted_response_kernels"].items()
                },
                "fingerprint_scope": "Complete rational coefficients rebuilt before exact monomial hashing; no finite-momentum truncation",
            }
        ),
        "whole_half_line_forcing_envelopes": serialize(source.envelopes()),
        "global_scalar_mean_and_observable_bounds": serialize(payload(bounds.data())),
        "independent_variable_and_asymptotic_controls": serialize(
            {
                "same_variable_change": identities.data()["same_mean_variable_change"],
                **payload(identities.asymptotic_controls()),
                "no_physical_divergence_inferred_from_failed_uncompensated_absolute_estimate": True,
            }
        ),
        "joint_seven_mode_relative_response": serialize(joint.data()),
        "explicit_amplitude_calibrations": serialize(
            {
                "scalar_only": calibration.calibrated(),
                "all_three_blocks": joint.calibrated(),
                "zero_addition": joint.calibrated(0, 0, 0),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "Actual weighted coupled scalar evolution and an exact compensated analysis variable yield integrable state-difference forcing on both infinite tails without changing the physical mean solution. The original trace is restored, and the actual physical matter observable has a uniform absolute t^-4 bound. The same reference with independent scalar, TT and Proca additions has globally controlled leading relative means and a complete strict-bounce exact-frame representative. Common absolute sources, co-evolved quantum geometry, nonlinear/higher-loop control, corrected cones and actual V/G/B matching remain open.",
        "not_established": [
            "A common absolute renormalized coupled source, reassignment of the older reference cancellation, or full interacting covariant/BRST Ward identity",
            "An exact finite-amplitude semiclassical solution, quantum state co-evolved on the representative, nonlinear/higher-loop or omitted-operator residual bounds",
            "Uniform fractional matter-density control, a conserved scalar-fluid identification, finite total R3 energy, or bounds on every Gaussian configuration",
            "Uniform corrected principal-cone/self-energy control, heavy thresholds/cutoff, actual vacuum/finite-gravity/common-parent V/G/B matching, UV completion or universal no-go",
            "Original P8(b) classification or original P8 closure",
        ],
        "verification_boundary": "Native rational phase/CCR and exact same-mean variable identities; every complete compensated kernel reconstructed; coefficientwise whole-half-line integrals and pointwise observable bounds; independent source parity, asymptotic cancellation and Cartesian scalar-TT contraction controls; exact amplitude and joint bounce inequalities. State, norm-ODE, block-source and geometric proofs are written/source-pinned, not proof-assistant formalized or peer reviewed. All ancestors, own source/proof/test bytes and report fields are read-only replayed using unmodified scientific SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The global seven-mode relative mean report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.107.GLOBAL_SEVEN_MODE_RELATIVE_MEAN_CONTROL replay passed; absolute source, higher orders, V/G/B and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
