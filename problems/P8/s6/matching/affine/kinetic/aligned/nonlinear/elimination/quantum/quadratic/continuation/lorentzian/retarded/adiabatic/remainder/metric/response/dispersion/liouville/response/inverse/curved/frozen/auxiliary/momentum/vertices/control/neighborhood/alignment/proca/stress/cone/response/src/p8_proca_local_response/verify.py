"""Read-only new finite local ordinary Proca response and fixed local-profile part."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as original
from p8_proca_scalar_cone import verify as parent

from . import audit, bounds, chart, local

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"local-ordinary-Proca-response.json"
PARENT_SHA="4b6ca08197fa6084c0e541ad4aa1dc4f786ce4673a212dca4f4238d1617ef407"
sha,serialize=original.sha,original.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_proca_local_response/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen central classical scalar cone report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_83_fully_rebuilt":PARENT_SHA,
            "classical_central_cone_not_promoted_to_a_quantum_cone":True}


@cache
def build_report():
    prior=prior_checks()
    identities,gates=audit.residuals(),audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A new local ordinary Proca response gate failed")
    exact=original.certify_residuals(identities)
    return {"schema":1,"claim":"P8-S6.84.LOCAL_ORDINARY_PROCA_RESPONSE","date":"2026-09-08",
            "status":"NEW_CONSTANT_MASS_FINITE_LOCAL_METRIC_RESPONSE_AND_FIXED_LOCAL_PROFILE_PART_CONTROLLED; FULL_NONLOCAL_INVERSE_MATCHING_AND_ORIGINAL_P8_OPEN",
            "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
            "formulation":"FORMULATION.md",
            "written_proofs":["notes/matching.md","notes/readouts.md","notes/chart.md",
                              "notes/bounds.md","notes/scope.md"],
            "exact_residuals":exact,"named_exact_check_count":len(identities),
            "checked_scalar_entries":sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in identities.values()),
            "proof_checks":gates,
            "new_finite_covariant_heat_quadratic_densities":serialize({order:local.density(order) for order in range(3)}),
            "actual_physical_coordinate_Euler_currents":serialize({
                output:{order:local.operator(output,order) for order in range(3)} for output in ("N","Z")}),
            "actual_normalized_local_physical_stress_response":serialize({
                output:{order:local.physical_operator(output,order) for order in range(3)} for output in ("energy","pressure")}),
            "actual_nonlinear_clock_chart_derivatives_and_time_jets":serialize(chart.mapping()),
            "actual_metric_pullback_and_fixed_local_profile_densities":serialize({
                order:chart.density(order) for order in range(3)}),
            "actual_clock_chart_local_Euler_currents":serialize({
                output:{order:chart.operator(output,order) for order in range(3)} for output in ("N","V")}),
            "new_rank_one_fourth_derivative_coefficient_and_null_direction":serialize(chart.top()),
            "continuous_local_response_coefficient_envelopes":serialize(bounds.envelopes()),
            "same_scale_new_local_response_C4_to_C0_bounds":serialize(bounds.at_scale(bounds.SCALE)),
            "controls":audit.controls(),
            "verdict":"The actual constant-mass ordinary Proca pole and finite local matching are independently recomputed with all mass-source jets zero. Physical readout normalization, the second nonlinear metric-map contact and the already fixed local-profile part are retained. The constant vacuum-density term cancels its fixed-profile counterpart. Both physical and clock-chart coordinate Hessians satisfy the appropriate weighted adjoint identities. The new local fourth-derivative clock matrix is nonzero rank one, not the old invertible indefinite matrix. All four continuous C4-to-C0 component bounds are below 10^-790 at fixed m=1000,L=10^400.",
            "not_established":["The full state-dependent nonlocal Gaussian response and the rest of the fixed profile",
                "A no-loss coupled causal inverse, all-frequency spectrum or quantum stability/cone",
                "Higher-derivative truncation resummation as an exact particle spectrum",
                "Canonical profile interactions, mixed loops, omitted operators or heavy-threshold error",
                "Finite Wilson matching, common-parent V/G/B, UV completion/exclusion or original P8 closure"],
            "verification_boundary":"Exact new dimensional/heat-action comparison, physical and nonlinear-chart variations, weighted Hessian adjoints and continuous rational coefficient bounds with explicit controls. The full parent chain is natively rebuilt. Finite local homogeneous response only; not proof-assistant formalization or a full quantum solution."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The local ordinary Proca response report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.84.LOCAL_ORDINARY_PROCA_RESPONSE replay passed; full response/inverse, matching and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()

