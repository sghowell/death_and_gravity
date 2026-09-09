"""Read-only actual coupled-scalar mean-source and compact-response certificate."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as certificate
from p8_tensor_mean_source import verify as parent

from . import audit, bounds, calibration, geometry, model, observable, response, source

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "actual-coupled-scalar-mean-response.json"
PARENT_SHA = "376eb054220d5c3aa9147de2a55f88043452aca53cb54591e992110a50942b68"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_scalar_mean_source/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def rational_signature(value):
    if isinstance(value, bool) or not isinstance(value, (int, sp.Expr)):
        raise TypeError("Require an exact scalar rational kernel")
    value = sp.sympify(value)
    if value.has(
        sp.Float, sp.oo, -sp.oo, sp.nan, sp.zoo, sp.I
    ) or not value.free_symbols <= {model.u, model.k}:
        raise ValueError(
            "Kernel fingerprints require real exact declared-variable expressions"
        )
    numerator, denominator = sp.fraction(sp.cancel(value))

    def signature(expr):
        poly = sp.Poly(expr, model.u, model.k, domain=sp.QQ)
        rows = [
            [list(powers), int(co.p), int(co.q)] for powers, co in sorted(poly.terms())
        ]
        encoded = json.dumps(rows, separators=(",", ":"), ensure_ascii=True).encode(
            "ascii"
        )
        return {
            "monomials": len(rows),
            "total_degree": int(poly.total_degree()),
            "sorted_exact_rational_monomial_sha256": hashlib.sha256(
                encoded
            ).hexdigest(),
        }

    return {"numerator": signature(numerator), "denominator": signature(denominator)}


def matrix_signature(value):
    return {
        "rows": value.rows,
        "columns": value.cols,
        "row_major_entries": [rational_signature(entry) for entry in value],
    }


def payload(data):
    return {key: value for key, value in data.items() if key != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen actual tensor-source certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_105_fully_rebuilt": PARENT_SHA,
        "same_original_retuned_action_and_physical_metric": True,
        "same_reference_other_five_modes_and_counterterms": True,
        "only_a_named_natural_scalar_band_covariance_added": True,
    }


@cache
def build_report():
    prior = prior_checks()
    identities = audit.residuals()
    gates = {key: bool(value) for key, value in audit.gates().items()}
    if not all(value is True for value in gates.values()):
        raise ValueError("An actual coupled-scalar source proof gate failed")
    exact = certificate.certify_residuals(identities)
    scalar = source.data()
    obs = observable.data()
    b = bounds.response_bounds()
    return {
        "schema": 1,
        "claim": "P8-S6.106.ACTUAL_COUPLED_SCALAR_MEAN_RESPONSE",
        "date": "2026-09-09",
        "status": "COMPLETE_ACTUAL_COUPLED_SCALAR_QUADRATIC_BACKGROUND_VARIATION_FOUR_MEAN_SOURCES_POSITIVE_FINITE_STATE_DIFFERENCE_AND_PHYSICAL_MATTER_OBSERVABLE_WITH_EXACT_COMPACT_BOUNCE_INTERVAL_CONTROL; NOT_UNIFORM_INFINITE_TAIL_CONTROL_AN_ABSOLUTE_OR_EXACT_SEMICLASSICAL_SOLUTION_HIGHER_ORDER_CONTROL_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/action.md",
            "notes/phase.md",
            "notes/mean.md",
            "notes/observable.md",
            "notes/bounds.md",
            "notes/center.md",
            "notes/scope.md",
        ],
        "exact_residuals": exact,
        "named_exact_check_count": len(identities),
        "checked_scalar_entries": sum(
            value.rows * value.cols if isinstance(value, sp.MatrixBase) else 1
            for value in identities.values()
        ),
        "proof_checks": gates,
        "actual_scalar_quadratic_background_action": serialize(
            {
                "literal_lapse_coefficient_jets": model.coefficients(),
                "generic_independent_mean_action": model.generic(),
                "natural_geometry_and_York": payload(geometry.data()),
            }
        ),
        "actual_natural_scalar_phase": serialize(
            {
                "canonical_map": payload(geometry.canonical()),
                "generator": matrix_signature(scalar["natural_phase_generator"]),
                "added_anchor_covariance": "eta*b(k)*I4 in natural phase, other five modes unchanged",
                "actual_covariance": "U_natural(u,0,k) eta*b(k)*I4 U_natural(s,0,k)^T",
                "band": "Same normalized smooth open shell 1<|k|^2<4; eta=hbar*lambda/kappa",
                "state_and_mean_all_finite_times_not_uniform_tail_bound": True,
            }
        ),
        "actual_leading_scalar_sources": serialize(
            {
                "definitions": {
                    "lapse": "A^-3 partial_N H2",
                    "scale": "(3A^3)^-1 partial_p H2",
                    "trace": "A^-3[-A partial_A H2/3+l partial_l H2]",
                    "matter_field": "A^-3 partial_l H2",
                },
                "complete_source_Hessian_fingerprints": {
                    key: matrix_signature(value)
                    for key, value in scalar["source_Hessians"].items()
                },
                "fingerprint_scope": "Exact numerator/denominator monomial hashes identify every complete rational kernel; source code reconstructs all coefficients before hashing",
            }
        ),
        "actual_leading_scalar_mean_response": serialize(payload(response.data())),
        "actual_physical_matter_observable": serialize(
            {
                "linear_inhomogeneous_lapse": obs["linear_inhomogeneous_lapse"],
                "additional_mean_density_and_pressure": obs[
                    "additional_mean_density_and_pressure"
                ],
                "intrinsic_density_Hessian": matrix_signature(
                    obs["intrinsic_density_Hessian"]
                ),
                "intrinsic_pressure_Hessian": matrix_signature(
                    obs["intrinsic_pressure_Hessian"]
                ),
                "independent_center_chart": payload(source.density_center()),
                "independent_full_center_observable": payload(observable.center()),
                "mean_lapse_counted_once": True,
            }
        ),
        "independent_scalar_center_jets": serialize(payload(response.center())),
        "whole_interval_state_and_response_bounds": serialize(
            {
                **payload(b),
                "actual_generator_row_sums": bounds.data()["generator_row_sums"],
                "actual_source_half_entry_sums": bounds.data()[
                    "source_half_entry_sums"
                ],
                "explicit_calibration": calibration.calibrated(),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "Independent homogeneous-background variation of the original scalar quadratic action yields the actual lapse, scale, trace and matter-field sources. Both canonical shifts and the time boundary reproduce the existing coupled scalar dynamics exactly. A named positive smooth scalar covariance addition preserves the reference's generalized Hadamard condition and its other five modes. Full center constraint/physical-density polynomials and exact whole-strip bounds control the leading mean response and actual matter observable. A positive local bouncing metric representative is obtained, but uniform scalar tails, an absolute or self-consistent quantum solution, higher-order errors and V/G/B matching remain open.",
        "not_established": [
            "Uniform infinite-tail response bounds, a globally complete scalar-perturbed representative, or a quantum state co-evolved on a new mean metric",
            "A common absolute coupled source, full covariant/BRST quantum Ward identity, exact finite-amplitude semiclassical solution, or nonlinear/higher-loop remainder control",
            "A conserved homogeneous fluid identification for the coordinate-averaged matter jet, finite total R3 energy, or positive lapse for every unbounded Gaussian field configuration",
            "A heavy spectrum/cutoff, omitted-operator bounds, actual vacuum/finite-gravity/common-parent V/G/B matching, UV completion, or a universal exclusion from one quadratic sign",
            "Original P8(b) classification or original P8 closure",
        ],
        "verification_boundary": "Native actual retuned coefficient reconstruction, independent canonical-volume/York/lapse elimination, mean variation before background restriction, complete time-dependent scalar phase bridge, source transport and full center observable/constraint comparisons. Whole-strip rational bounds reject unsupported coefficient structures. State/ODE and local geometric proofs are source-pinned but not proof-assistant formalized. Complete kernel fingerprints are recomputed from exact rational coefficients, not substitutes for kernel derivation. Every ancestor, own source/proof/test byte and report field is independently replayed with unmodified scientific SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The actual coupled-scalar mean-source report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.106.ACTUAL_COUPLED_SCALAR_MEAN_RESPONSE replay passed; uniform scalar tails, absolute source, V/G/B and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
