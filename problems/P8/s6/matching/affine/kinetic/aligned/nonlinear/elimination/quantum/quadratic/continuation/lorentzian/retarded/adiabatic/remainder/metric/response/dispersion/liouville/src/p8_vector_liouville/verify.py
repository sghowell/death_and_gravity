"""Read-only exact vector acoustic-time reduction and prepared potential bound."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_vector_metric_dispersion import verify as parent

from . import clock, proofs, reduction, variation

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"acoustic-vector-reduction.json"
PARENT_SHA = "decbb30e4afef9c5f045678b664162c9af3f908c21cffd8f1a0812287e3c36f3"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_vector_liouville/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen metric dispersion and chart report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_69_fully_rebuilt": PARENT_SHA, "physical_metric_action_state_and_subtraction_unchanged": True}


def controls():
    calls = []
    for value in (True, False, sp.true, sp.false, 0, 1, 1.0, sp.Integer(1), "transverse", "", None):
        calls.append(lambda value=value: reduction.data(value))
    for value in (True, False, sp.true, sp.false, 1.0, sp.Float(1), "1", sp.I, sp.oo, -sp.oo, sp.zoo, sp.nan, -1, None):
        calls.extend((lambda value=value: clock.bounds(value, 1000),
                      lambda value=value: clock.bounds(1, value),
                      lambda value=value: variation.bound(value, 1000),
                      lambda value=value: variation.bound(1, value)))
    calls.extend((lambda: clock.bounds(1, 0), lambda: variation.bound(1, 0)))
    for value in (0, 1, "yes", None, sp.true):
        calls.append(lambda value=value: clock.exact_rational(1, value))
    for value in (True, False, sp.true, sp.false, 1.0, sp.Float(1), "0", None,
                  variation.n[3], variation.zeta[3], variation.n[0]*variation.zeta[0],
                  1, sp.Symbol("foreign")*variation.n[0], sp.oo, sp.nan):
        calls.append(lambda value=value: variation.coefficient_bound(value))
    calls.extend(lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan))
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An unsupported acoustic reduction input was accepted")
    return {"rejected_inputs": rejected, "native_domain_validation_before_cached_implementations": True,
            "both_third_source_derivatives_excluded_from_C2": True,
            "zero_polynomial_residual_accepted_without_accepting_nonzero_constants": True,
            "physical_readout_acoustic_rate_retained": True,
            "prepared_coordinate_shift_and_initial_history_retained": True,
            "momentum_zero_is_a_potential_limit_not_an_invertible_mode_map": True}


@cache
def build_report():
    prior = prior_checks()
    identities, gates = proofs.residuals(), proofs.checks()
    if not all(value is True for value in gates.values()):
        raise ValueError("A vector acoustic-time audit gate failed")
    certified = affine.certify_residuals(identities)
    return {"schema": 1, "claim": "P8-S6.70.ACOUSTIC_VECTOR_REDUCTION", "date": "2026-09-08",
            "status": "EXACT_ACOUSTIC_MODE_REDUCTION_AND_PREPARED_POTENTIAL_BOUND; ORIGINAL_P8_OPEN",
            "prior_sha256": prior, "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
            "formulation": "FORMULATION.md", "written_proofs": ["notes/reduction.md", "notes/variation.md", "notes/boundary.md"],
            "exact_residuals": certified, "named_exact_check_count": len(identities),
            "checked_scalar_entries": sum(v.rows*v.cols if isinstance(v, sp.MatrixBase) else 1 for v in identities.values()),
            "proof_checks": gates,
            "positive_physical_mode_reduction": serialize({sector: reduction.data(sector) for sector in ("T", "L")}),
            "longitudinal_scalar_and_momentum_potentials": serialize(reduction.longitudinal_potential()),
            "original_clock_coefficients": serialize(clock.data()),
            "prepared_source_chart": serialize(variation.source()),
            "prepared_potential_variation": serialize(variation.data()),
            "continuous_prepared_variation_majorants": serialize(variation.estimates()),
            "fixed_scale_examples": serialize({str(k): {"background": clock.bounds(k, 1000),
                                                        "prepared_variation": variation.bound(k, 1000)}
                                                for k in (0, 1000, 10000)}),
            "controls": controls(),
            "verdict": "Both physical mode sectors admit exact positive Liouville clocks preserving Wronskians and physical currents. The longitudinal pump yields a scalar-shaped principal potential plus an explicit finite-momentum correction. On I the latter lies between zero and 88*m^2/(k^2+m^2); its prepared fixed-acoustic-time first variation has C2-to-C0 norm below 18000*m^2/(k^2+m^2), retaining the retarded coordinate shift. The regular pump/rate source chart reproduces the previous rank-one UV pole direction.",
            "not_established": ["An exact replacement scalar theory, off-clock ordinary-Proca state/cone transfer or an invertible canonical k=0 map",
                                "The full renormalized stress decomposition or a no-loss coupled causal inverse",
                                "Spatial response, independent initial states, quantum stability/cones, interactions, other loops, cutoff, finite Wilson matching, V/G/B or original P8 closure"],
            "verification_boundary": "Exact canonical, time-change, potential, variation and source-chart identities with written continuous momentum and source bounds. Potential-level, not stress-level. Not proof-assistant formalization."}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The acoustic vector-reduction report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.70.ACOUSTIC_VECTOR_REDUCTION replay passed; coupled response and original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
