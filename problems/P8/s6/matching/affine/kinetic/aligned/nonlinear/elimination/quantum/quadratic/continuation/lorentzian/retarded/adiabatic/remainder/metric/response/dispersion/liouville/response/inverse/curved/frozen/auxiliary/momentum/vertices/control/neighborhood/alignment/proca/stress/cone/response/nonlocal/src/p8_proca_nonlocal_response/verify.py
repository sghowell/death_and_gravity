"""Read-only complete prepared ordinary Proca response with its fixed profile."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as original
from p8_proca_local_response import bounds as local_bounds
from p8_proca_local_response import verify as parent

from . import audit, envelopes, evolution, source, tadpole, tail, tangent

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"prepared-ordinary-Proca-response.json"
PARENT_SHA="240b1a0550969a51bd1c0ea4043360566a1b1cd10408ea30b2f89d560dc1aea2"
sha,serialize=original.sha,original.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_proca_nonlocal_response/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen new local ordinary Proca response report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_84_fully_rebuilt":PARENT_SHA,
            "new_nonlocal_response_recomputed_not_old_metric_mass_response_transferred":True}


@cache
def build_report():
    prior=prior_checks()
    identities,gates=audit.residuals(),audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A new prepared ordinary Proca response gate failed")
    exact=original.certify_residuals(identities)
    return {"schema":1,"claim":"P8-S6.85.PREPARED_ORDINARY_PROCA_RESPONSE","date":"2026-09-08",
            "status":"NEW_COMPLETE_PREPARED_HOMOGENEOUS_GAUSSIAN_PLUS_FIXED_PROFILE_RESPONSE_CONTROLLED_IN_BOTH_CHARTS; COUPLED_INVERSE_MATCHING_AND_ORIGINAL_P8_OPEN",
            "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
            "formulation":"FORMULATION.md",
            "written_proofs":["notes/vertices.md","notes/reference.md","notes/transport.md",
                              "notes/profile-chart.md","notes/limit.md"],
            "exact_residuals":exact,"named_exact_check_count":len(identities),
            "checked_scalar_entries":sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in identities.values()),
            "proof_checks":gates,
            "actual_new_mode_source_vertices":serialize({sector:source.data(sector) for sector in ("T","L")}),
            "actual_new_normalized_physical_readout_vertices":serialize({
                sector:{output:source.readout(sector,output) for output in ("energy","pressure")} for sector in ("T","L")}),
            "actual_new_varied_adiabatic_readouts":serialize({
                sector:{output:{order:source.adiabatic(sector,output,order) for order in range(3)}
                        for output in ("energy","pressure")} for sector in ("T","L")}),
            "new_reference_low_residuals_and_full_tenth_source_coefficients":serialize({
                sector:{key:tangent.reference(sector)[key] for key in
                        ("varied_residual_low_coefficients","tenth_source_derivative_coefficients","tenth_source_derivative_fixtures")}
                for sector in ("T","L")}),
            "new_continuous_varied_reference_envelopes":serialize({sector:envelopes.reference(sector) for sector in ("T","L")}),
            "new_integrable_reference_readout_tails":serialize({
                sector:{output:tail.reference_tail(sector,output) for output in ("energy","pressure")} for sector in ("T","L")}),
            "new_transport_constants_with_original_nonzero_initial_mixing":serialize(evolution.constants()),
            "new_complete_physical_vector_response_bounds":serialize(evolution.bound(local_bounds.SCALE)),
            "actual_already_fixed_profile_vertices_and_total_chart_contact":serialize(tadpole.physical_vertices()),
            "continuous_actual_clock_chart_C10_lift":serialize(tadpole.chart_lift()),
            "new_complete_background_cancelled_response_bounds_in_both_charts":serialize(tadpole.bound(local_bounds.SCALE)),
            "controls":audit.controls(),
            "verdict":"The new constant-mass prepared homogeneous response is rebuilt from literal physical vertices, exact covariance transport, newly varied eighth-order reference and integrable fourth-order subtraction. Both the exact-mode and all six varied-subtraction Ward identities pass. Original nonzero initial mixing is retained. Adding the independently matched new local operator and the full already-fixed S6.82 profile gives a physical C10-to-C0 bound below 10^-790 at fixed m=1000,L=10^400; the actual clock-chart bound is below 10^-779. The second nonlinear metric-map contact cancels only for the total background-cancelled action. The common-dimensional limit and retardedness hold on the declared prepared homogeneous source domain.",
            "not_established":["A no-loss coupled causal inverse, feedback contraction or all-frequency quantum stability/cone",
                "Arbitrary spatial or independently varied initial-state response",
                "A nonlinear quantum bounce or canonical interaction bounds for the fixed profile",
                "Mixed loops, omitted operators, heavy thresholds or an interacting Wilsonian cutoff",
                "Finite Wilson matching, common-parent V/G/B, UV completion/exclusion or original P8 closure"],
            "verification_boundary":"Exact new physical-mode/contact/Ward/subtraction identities, continuous source coefficient and full radial transport bounds, and written common-limit, causal and chart arguments. Every parent is natively rebuilt. Prepared homogeneous retained Gaussian response with C10-to-C0 derivative loss, not proof-assistant formalization or a full quantum solution."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The prepared ordinary Proca response report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.85.PREPARED_ORDINARY_PROCA_RESPONSE replay passed; coupled inverse, matching and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()

