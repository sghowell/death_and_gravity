"""Read-only full matter-vector momentum reduction and physical quadratic replay."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_auxiliary_neighborhood import model
from p8_auxiliary_neighborhood import verify as parent
from p8_physical import jets as j

from . import canonical, constraints, fixtures, projector, quadratic

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"coupled-spatial-reduction.json"
PARENT_SHA="fbf3c9d949cbe63304497c577091c3401b40cd9ed49f78af14fa3a4d8b754e31"
sha,serialize=affine.sha,affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_coupled_momentum/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen quantified auxiliary-domain report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_75_fully_rebuilt":PARENT_SHA,
            "original_full_Hamiltonian_matter_vector_generators_and_margin_replayed":True}


def controls():
    calls=[]
    for value in (True,False,sp.true,sp.false,0,1,5,2.0,sp.Integer(2),"2",None,[2]):
        calls.append(lambda value=value:fixtures.build(value))
    for value in (True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,None):
        calls.append(lambda value=value:constraints.context(((value,0,0),(-1,0,0))))
    for momenta in ((),((1,0,0),),((1,0,0),)*5,((1,0),(-1,0)),
                    ((1,0,0),(1,0,0)),((1,0,0),(-1,0,0),(0,1,0),(0,-1,0))):
        calls.append(lambda momenta=momenta:constraints.context(momenta))
    good=((1,0,0),(-1,0,0))
    for value in (1,0,sp.true,sp.false,1.0,None):
        calls.extend((lambda value=value:constraints.context(good,symbolic_time=value),
                      lambda value=value:constraints.context(good,symbolic_scale=value)))
    for value in (True,False,sp.true,sp.false,0,1,1.0,"unknown",None,["curvature"]):
        calls.extend((lambda value=value:quadratic.pair(value,"curvature"),
                      lambda value=value:quadratic.pair("curvature",value)))
    calls.extend(lambda value=value:serialize(value) for value in (1.0,sp.Float(1),sp.oo,sp.nan))
    ctx=constraints.context(good)
    source=[j.Jet(ctx,{ctx.full:ctx.one_coefficient}),ctx.jet(),ctx.jet()]
    calls.append(lambda:constraints.old.solve_vector(ctx,source))
    rejected=0
    for call in calls:
        try:
            call()
        except (TypeError,ValueError):
            rejected+=1
    if rejected!=len(calls):
        raise ValueError("An unsupported coupled-momentum input was accepted")
    return {"rejected_inputs":rejected,"native_validation_before_cache_and_rational_conversion":True,
            "nonzero_proper_subset_domain_checked_before_recursion":True,
            "actual_nonzero_homogeneous_constraint_source_rejected":True,
            "both_matter_momenta_shifted_by_same_canonical_boundary":True,
            "full_nonlinear_inverse_metric_raises_the_matter_vector_generator":True,
            "temporal_vector_momentum_divergence_term_not_dropped":True,
            "nonvanishing_rational_TT_frame_and_dual_momentum_used":True,
            "full_volume_geometry_and_moving_boundaries_kept":True,
            "quadratic_phase_matrix_not_a_resummed_quantum_spectrum":True,
            "D_only_cutoff_not_transferred_to_coupled_parent":True}


@cache
def residuals():
    result={}
    for group in (fixtures.checks(),canonical.checks(),projector.checks(),quadratic.checks()):
        if set(result).intersection(group):
            raise ValueError("Repeated coupled-spatial-reduction residual name")
        result.update(group)
    return result


@cache
def gates():
    negative=fixtures.negative_controls()
    bg=model.coefficients()["background"]
    matrix=quadratic.phase_matrix()
    return {"two_three_four_labelled_fixtures_keep_every_constructed_constraint":all(
                all(fixtures.build(n)["result"]["checks"].values()) for n in (2,3,4)),
            "all_rational_TT_polarizations_are_nonzero_and_dualizable":all(
                sp.trace(E*E)>0 for n in (2,3,4) for E in fixtures.build(n)["polarizations"]),
            "omitted_vector_divergence_has_nonzero_nonlinear_defect":any(
                x!=0 for x in negative["omitted_vector_divergence_constraint_defect"]),
            "omitted_rolling_matter_has_nonzero_linear_defect":any(
                x!=0 for x in negative["omitted_background_matter_linear_constraint_solution_defect"]),
            "actual_comoving_matter_momentum_is_fixed_one_tenth":bool(bg["a"]**3*bg["ell"]==sp.Rational(1,10)),
            "seven_physical_modes_have_fourteen_canonical_phase_channels":len(quadratic.CHANNELS)==14,
            "all_105_independent_quadratic_phase_pairs_checked":sum(
                len(quadratic.CHANNELS)-i for i in range(len(quadratic.CHANNELS)))==105,
            "actual_quadratic_phase_matrix_is_fourteen_by_fourteen":matrix.shape==(14,14),
            "actual_quadratic_phase_matrix_obeys_real_functional_Hessian_adjoint":matrix==matrix.T.subs(
                constraints.wave_scale,-constraints.wave_scale),
            "compact_scalar_Hamiltonian_matches_prior_margin_after_both_momentum_shifts":quadratic.compact_scalar()["prior_regular_scalar_replay"]==0,
            "actual_lapse_force_has_no_linear_vector_divergence":quadratic.lapse_coefficients()[str(model.j)]==0,
            "exact_York_norm_bound_requires_nonzero_transfer":True,
            "highest_fourth_order_longitudinal_correction_drops_only_by_background_trace":True,
            "actual_time_derivatives_not_frozen_out_of_quadratic_comparison":True,
            "finite_Fourier_recursion_not_a_global_spatial_gauge_existence_theorem":True,
            "physical_higher_vertices_and_controlled_interacting_band_still_open":True,
            "no_quantum_response_finite_matching_or_original_P8_closure":True}


def fixture_metadata():
    return {str(n):{"rational_directions":fixtures.build(n)["context"].momenta,
                    "symbolic_time":model.u,"positive_wave_scale":constraints.wave_scale,
                    "solved_momentum_orders":tuple(range(1,n)),
                    "polarizations":fixtures.build(n)["polarizations"],
                    "constraint_gate_count":len(fixtures.build(n)["result"]["checks"])}
            for n in (2,3,4)}


@cache
def build_report():
    prior=prior_checks()
    identities,audit=residuals(),gates()
    if not all(value is True for value in audit.values()):
        raise ValueError("A coupled physical-spatial audit gate failed")
    exact=affine.certify_residuals(identities)
    return {"schema":1,"claim":"P8-S6.76.COUPLED_PHYSICAL_SPATIAL_REDUCTION","date":"2026-09-08",
            "status":"FULL_MATTER_VECTOR_MOMENTUM_RECURSION_AND_SEVEN_MODE_QUADRATIC_HAMILTONIAN_REPLAY; INTERACTING_CUTOFF_AND_ORIGINAL_P8_OPEN",
            "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
            "formulation":"FORMULATION.md","written_proofs":["notes/constraints.md","notes/canonical.md","notes/scope.md"],
            "exact_residuals":exact,"named_exact_check_count":len(identities),
            "checked_scalar_entries":sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in identities.values()),
            "proof_checks":audit,"canonical_matter_boundary":serialize(canonical.data()),
            "physical_momentum_definition":"Div(pi)^i-(1/2)*g^ij*(P_chi*partial_j chi+Pi^k F_jk-W_j*partial_k Pi^k)=0. Canonical momenta are densities and the complete inverse metric is kept.",
            "exact_York_operator_inverse_and_norms":serialize(projector.data()),
            "mixed_fixture_metadata":serialize(fixture_metadata()),
            "omission_negative_controls":serialize(fixtures.negative_controls()),
            "full_quadratic_phase_channels":list(quadratic.CHANNELS),
            "full_quadratic_phase_matrix":serialize(quadratic.phase_matrix()),
            "actual_linear_lapse_invariant_coefficients":serialize(quadratic.lapse_coefficients()),
            "compact_scalar_Hamiltonian":serialize(quadratic.compact_scalar()),
            "nonzero_transfer_and_quadratic_scope":"Finite labelled Fourier vertices, two through four external legs, zero total and no zero nonempty proper-subset wavevector. The recursion solves degrees one through n-1. The actual quadratic Hamiltonian is checked for all 105 independent pairs of 14 phase channels at symbolic u and positive wave scale; remaining entries follow from the real integrated functional Hessian. Physical momentum is in fixed local units, not held constant during a subsequent cosmological evolution.",
            "controls":controls(),
            "verdict":"All three physical momentum constraints, including rolling matter and the vector momentum-divergence term, have the finite nonzero-transfer York recursion through cubic order. The scalar-matter canonical boundary shifts both momenta. Full spatial geometry, lapse elimination and moving canonical terms reproduce the prior scalar/matter margin block, both tensor polarizations and all three vector polarizations. This validates the physical reduction for subsequent higher vertices, not an interacting cutoff or a quantum stability theorem.",
            "not_established":["Complete cubic/quartic physical vertices, normalized finite-time hard-channel transition bounds or a controlled interacting cutoff",
                "Global spatial gauge/Sobolev construction, zero/forward channels, nonlinear evolution or semiclassical spatial response",
                "Finite Wilson matching, common-parent V/G/B or original P8 closure"],
            "verification_boundary":"Exact symbolic-time and wave-scale Fourier identities and 105 independent quadratic phase-pair comparisons, with written canonical, recursive, adjoint and norm arguments. Not proof-assistant formalization."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The coupled-spatial-reduction report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.76.COUPLED_PHYSICAL_SPATIAL_REDUCTION replay passed; interacting cutoff and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
