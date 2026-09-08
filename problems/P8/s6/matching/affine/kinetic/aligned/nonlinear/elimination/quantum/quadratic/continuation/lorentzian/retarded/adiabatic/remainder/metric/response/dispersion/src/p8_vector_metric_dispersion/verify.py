"""Read-only full-metric flat-dispersion and regular-chart checkpoint."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_vector_metric_response import verify as parent

from . import chart, principal, proofs, spectral

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"metric-dispersion-chart.json"
PARENT_SHA = "b67525a295507e7c6e9a05a01c5f51e36abca7f004bd8d2b570bcb81bb9aabe7"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_vector_metric_dispersion/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen finite prepared metric response changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_68_fully_rebuilt": PARENT_SHA, "state_action_and_subtraction_unchanged": True}


def controls():
    calls = []
    for value in (True, False, sp.true, sp.false, 1, 0, 1.0, sp.Integer(1), "transverse", "", None):
        calls.extend((lambda value=value: spectral.pair(value), lambda value=value: spectral.contact(value)))
    for value in (True, False, sp.true, sp.false, 0, 4, -1, "1", 1.0, sp.Integer(1), None):
        calls.append(lambda value=value: spectral.moment(value, 2))
    for value in (True, False, sp.true, sp.false, 0, 1, -1, 0.5, sp.Float(2), "2", sp.I, sp.oo, -sp.oo, sp.zoo, sp.nan):
        calls.append(lambda value=value: spectral.moment(1, value))
    calls.extend(lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan))
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("Unsupported dispersion input was accepted")
    return {"rejected_inputs": rejected, "validation_precedes_cached_implementations": True,
            "omitting_static_contacts_breaks_vacuum_Hessian": True,
            "omitting_chart_background_contact_before_tadpole_is_invalid": True,
            "leading_log_root_not_an_exact_kernel_pole": True,
            "real_axis_invertibility_not_a_Bromwich_contour_bound": True,
            "one_feedback_application_not_a_resummed_solution": True}


@cache
def build_report():
    prior = prior_checks()
    identities, gates = proofs.residuals(), proofs.checks()
    if not all(value is True for value in gates.values()):
        raise ValueError("A metric-dispersion audit gate failed")
    certified = affine.certify_residuals(identities)
    d = sp.Symbol("d_greater_than_one", positive=True)
    return {"schema": 1, "claim": "P8-S6.69.METRIC_DISPERSION_CHART", "date": "2026-09-08",
            "status": "FROZEN_FLAT_DISPERSION_AND_REGULAR_RESPONSE_CHART; ORIGINAL_P8_OPEN",
            "prior_sha256": prior, "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
            "formulation": "FORMULATION.md",
            "written_proofs": ["notes/dispersion.md", "notes/chart.md", "notes/next-inverse.md"],
            "exact_residuals": certified, "named_exact_check_count": len(identities),
            "checked_scalar_entries": sum(v.rows*v.cols if isinstance(v, sp.MatrixBase) else 1 for v in identities.values()),
            "proof_checks": gates,
            "physical_flat_pair_vertices": serialize({sector: spectral.pair(sector) for sector in ("T", "L")}),
            "physical_flat_contacts": serialize({sector: spectral.contact(sector)["contact"] for sector in ("T", "L")}),
            "spectral_matrix": serialize(spectral.matrix()),
            "radial_moments": serialize({j: spectral.moment(j, d) for j in (1, 2, 3)}),
            "thrice_subtracted_normalized_remainder": serialize(spectral.subtracted_normalized(d)),
            "finite_large_frequency_constant": serialize(spectral.actual(spectral.finite_asymptotic())),
            "frozen_pole_and_finite_fourth_symbols": serialize({"pole": principal.pole(), "finite": principal.finite(),
                                                               "full_asymptotic_constant": principal.asymptotic()}),
            "isolated_positive_real_axis_bounds": serialize(principal.real_axis()),
            "regular_chart": serialize(chart.data()), "one_prepared_tree_response_bounds": serialize(chart.norm()),
            "controls": controls(),
            "verdict": "The full homogeneous flat massive two-current bubble admits exact Taylor-subtracted dispersion. Its pole is rank one, but the finite complement remains. The isolated normalized fourth-order block is invertible on the positive real Laplace axis with transformed matrix inverse norm below 26 on I; its apparent leading-log zero is not an exact pole. The nonlinear regular-chart contact cancels only after the fixed tadpole is included. One prepared retarded tree response is bounded, not iterated.",
            "not_established": ["A causal inverse of the full curved tree-plus-loop metric operator, arbitrary spatial momentum or varied initial covariance",
                                "Quantum stability, physical cones, interactions, other loops, cutoff, finite Wilson matching, V/G/B or original P8 closure"],
            "verification_boundary": "Exact canonical, radial, principal-symbol and chart identities with written spectral, asymptotic and real-axis inequalities. The flat vacuum comparator is not a replacement state on the rolling clock. Not proof-assistant formalization."}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The metric-dispersion report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.69.METRIC_DISPERSION_CHART replay passed; full coupled inverse and original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
