"""Read-only complete coupled higher-vertex and continuous raw-bound replay."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_coupled_momentum import verify as parent

from . import anchors, audit, coefficients, majorant, physical, taylor

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"coupled-higher-vertices.json"
PARENT_SHA="5d5ceba34f20d2b7e2db2ed229f45746bb6ea208238f69ebb05736b44170bceb"
sha,serialize=affine.sha,affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_coupled_vertices/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen coupled spatial-reduction report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_76_fully_rebuilt":PARENT_SHA,
            "original_full_classical_parent_and_seven_mode_quadratic_matrix_replayed":True}


def fixture_reports():
    result={}
    for n,point in audit.POINTS:
        data=physical.fixture(n,point)
        result[f"{n}_legs_at_{point}"]={
            "time":point,"directions":data["context"].momenta,
            "kernel":data["full_labelled_kernel"],
            "spatial_constraint_checks":data["momentum"]["checks"],
            "physical_lapse_and_stationary_checks":data["reduction"]["checks"]}
    return result


@cache
def build_report():
    prior=prior_checks()
    identities,gates=audit.residuals(),audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A coupled higher-vertex gate failed")
    exact=affine.certify_residuals(identities)
    boundary=coefficients.boundary_jets()
    rows={",".join(map(str,key)):value for key,value in coefficients.jets().items()}
    return {"schema":1,"claim":"P8-S6.77.COUPLED_PHYSICAL_HIGHER_VERTICES","date":"2026-09-08",
            "status":"ACTUAL_SEVEN_MODE_CLASSICAL_PHYSICAL_H3_H4_AND_CONTINUOUS_RAW_VERTEX_BOUNDS; NORMALIZED_INTERACTION_CUTOFF_AND_ORIGINAL_P8_OPEN",
            "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
            "formulation":"FORMULATION.md","written_proofs":["notes/reduction.md","notes/bounds.md","notes/scope.md"],
            "exact_residuals":exact,"named_exact_check_count":len(identities),
            "checked_scalar_entries":sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in identities.values()),
            "proof_checks":gates,"actual_lapse_monomial_derivatives":serialize(rows),
            "original_boundary_and_Q_lapse_jets":serialize({
                name:taylor.derivatives(boundary[name]) for name in ("I","Iphi","Q")}),
            "actual_time_coefficient_bounds":serialize(majorant.coefficient_bounds()),
            "physical_vertex_fixture_kernels":serialize(fixture_reports()),
            "independent_electric_quartic":serialize(anchors.pure_electric()),
            "reality_and_permutation_fixtures":serialize({
                f"{n}_legs":{name:anchors.transformed_fixture(n,**kwargs) for name,kwargs in (
                    ("reality",{"reverse":True}),
                    ("permutation",{"permutation":tuple(reversed(range(n)))}))}
                for n in (3,4)}),
            "raw_physical_majorant_example":serialize(majorant.serialize_bound(majorant.bound(1,4,1))),
            "hard_wave_raw_majorant_example":serialize(majorant.serialize_bound(
                majorant.bound(1,10**9,sp.Rational(1,10**8)))),
            "vertex_domain_and_normalization":"I=[-1/2,1/2], fixed local units, two through four labelled legs with zero total and no zero nonempty proper-subset transfer. Raw canonical phase components, full factorial multiplicities and all matter/vector density factors. The continuous raw bound is conditional on specified component/derivative/inverse-transfer upper bounds, not on a list of sampled momenta or normalized oscillator columns.",
            "actual_stationary_formulas":"H=A(N)+L(N,Y)+Q(N,Y), A1=0,A2=-2Je. n1=-L1/A2,F2=A3*n1^2/2+L2*n1+Q1,n2=-F2/A2. H2=Q0-L1^2/(2A2); H3=A3*n1^3/6+L2*n1^2/2+Q1*n1; H4=A4*n1^4/24+L3*n1^3/6+Q2*n1^2/2-F2^2/(2A2). All physical geometry is substituted before extracting the labelled coefficient; highest fourth-order York momentum drops by background trace.",
            "controls":audit.controls(),
            "verdict":"The actual classical seven-mode spatial and auxiliary reduction defines H3/H4 on the finite nonzero-transfer Fourier domain. All 70 lapse derivatives retain the original time primitive and Q source; an all-time invariant identity transfers the full prior quadratic block. Complete mixed fixtures, an independent nonzero electric quartic and Fourier reality/permutation checks replay. A written positive recursion supplies continuous raw phase-component H3/H4 bounds. This is not a normalized transition norm, interacting cutoff, quantum interaction control or original P8 closure.",
            "not_established":["Normalized seven-mode finite-time transition or connected exchange bounds and controlled interacting cutoff",
                "Omitted loops/operators/thresholds, nonlinear semiclassical or global spatial existence, zero/forward channels",
                "Finite Wilson matching, common-parent V/G/B or original P8 closure"],
            "verification_boundary":"Exact symbolic lapse/invariant identities and rational Fourier fixtures, with written generic spatial/stationary reduction and continuous positive-majorant proofs. Not proof-assistant formalization."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The coupled higher-vertex report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.77.COUPLED_PHYSICAL_HIGHER_VERTICES replay passed; normalized cutoff and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
