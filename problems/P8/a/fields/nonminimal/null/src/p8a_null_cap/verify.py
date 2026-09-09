"""Read-only actual squeezed-state null-line obstruction certificate."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine.verify import serialize
from p8a_nonminimal import verify as prior

from . import audit, construction, independent, normalization

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "null-one-sided-state-cap-obstruction.json"
PARENT_SHA = "dfbb9ce85e818e50280aa838e43893f381f30f4573e879e9f820bc4b1e7b5292"
CLAIM = "P8-A.21.NULL_ONE_SIDED_STATE_CAP_OBSTRUCTION"
REPORT_KEYS = (
    "schema",
    "claim",
    "date",
    "status",
    "prior_sha256",
    "source_sha256",
    "formulation",
    "written_proofs",
    "exact_residuals",
    "named_exact_check_count",
    "checked_scalar_entries",
    "derived_constants",
    "proof_checks",
    "controls",
    "actual_state_and_normalization",
    "null_inequality_obstruction",
    "energy_and_ANEC",
    "verification_boundary",
    "not_established",
    "verdict",
)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def source_files():
    return (
        sorted(ROOT.glob("src/p8a_null_cap/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


@cache
def prior_checks():
    if sha(prior.REPORT) != PARENT_SHA:
        raise ValueError("The pinned physical nonminimal stress source changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {
        "A20_and_its_ancestry_fully_rebuilt": PARENT_SHA,
        "parent_arbitrary_jet_stress_used_for_actual_null_contraction": True,
        "parent_timelike_QSEI_applied_to_null_line": False,
    }


def exact_checks():
    rows = {**construction.checks(), **normalization.data()["checks"]}
    for name, v in rows.items():
        if sp.simplify(v) != 0:
            raise ValueError("An actual null-state identity failed: " + name)
    return dict.fromkeys(rows, "0")


def checked_constants():
    primary = {
        "profile_bounds": construction.profile_bounds(),
        "resource_bounds": normalization.resource_bounds(),
    }
    alternate = independent.replay()
    if serialize(primary) != serialize(alternate):
        raise ValueError("The independent Fraction null-state reconstruction differs")
    return {**serialize(primary), "independent_Fraction_reconstruction_agrees": True}


def build_report():
    inputs = prior_checks()
    rows = exact_checks()
    proofs = {k: bool(v) for k, v in audit.proof_checks().items()}
    if not all(proofs.values()):
        raise ValueError("A continuous null-state proof gate failed")
    result = {
        "schema": 1,
        "claim": CLAIM,
        "date": "2026-09-09",
        "status": "COMPLETE_EXPLICIT_ACTUAL_FINITE_ENERGY_HADAMARD_SQUEEZED_STATE_NULL_SEGMENT_OBSTRUCTION_EVEN_WITH_ONE_SIDED_WICK_UPPER_CAP_FOR_XI_ZERO_TO_QUARTER; FULL_LINE_ANEC_POSITIVE_FOR_THIS_FAMILY_GENERAL_NULL_THEOREM_AND_ORIGINAL_P8_OPEN",
        "prior_sha256": inputs,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": ["notes/state.md", "notes/profile.md", "notes/scope.md"],
        "exact_residuals": rows,
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": len(rows),
        "derived_constants": checked_constants(),
        "proof_checks": proofs,
        "controls": {
            "rejected_inputs": audit.rejected_inputs(),
            "one_sided_upper_cap_not_an_absolute_Wick_square_cap": True,
            "compact_negative_average_not_negative_complete_line_ANEC": True,
        },
        "actual_state_and_normalization": {
            "field": "one real massless scalar, physical improved stress, 0<=xi<=1/4 in four-dimensional Minkowski space",
            "null_line": "gamma(s)=(s,0,0,s), ell=(1,0,0,1), |s|<=1 in fixed affine sampling units",
            "physical_mass_shell": "omega=(p^2+q^2)/(2*q), kz=(p^2-q^2)/(2*q), q>0",
            "exact_measure": "dmu=dq*d^2p/(16*pi^3*q)",
            "fixed_frequency_data": "epsilon=1/1000,width=1/10000; even nonnegative C-infinity bump b of mass one",
            "F": "(4-i/epsilon)*b(q-epsilon)+(-2+i/epsilon)*b(q-2*epsilon)",
            "transverse_data": "fixed nonzero nonnegative C-infinity chi supported on 1<|p|<2; chi(p/R), R>=1",
            "actual_unit_mode": "h_R=q*F(q)*chi(p/R)/(R*sqrt(C)), C=Iq*J2/(16*pi^3)",
            "mode_integrals": "Ij=integral q^j*|F|^2 (Iq has j=1); J1=integral chi,J2=integral chi^2,Jp=integral p^2*chi^2",
            "actual_affine_solution": "u_R(gamma(s))=A_R*v(s), A_R=sqrt(hbar)*R*gamma0, gamma0=J1/(16*pi^3*sqrt(C))>0",
            "v": "bhat(s)*[(4-i/epsilon)*exp(-i*epsilon*s)+(-2+i/epsilon)*exp(-2*i*epsilon*s)]",
            "positive_pure_mode_squeeze": "n=1/3,m=-2/3; quadrature variances 1/6,3/2; vacuum on orthogonal complement",
            "Hadamard": "each finite R has compact smooth mass-shell data; the exact covariance difference is a smooth finite sum of mode outer products",
            "renormalized_smooth_difference_formed_before_null_restriction": True,
            "vacuum_distribution_itself_pulled_back_to_null_line": False,
        },
        "null_inequality_obstruction": {
            "full_physical_operator": "T_ll=(D_l Phi)^2-xi*D_l^2(Phi^2); improvement retained",
            "whole_segment_Wick_upper": "w_R/A_R^2<=-1633/2500<-3/5",
            "whole_segment_null_upper_all_xi": "T_ll/A_R^2<=-4747/15000<-3/10 for all 0<=xi<=1/4",
            "whole_segment_conformal_null_upper": "T_ll/A_R^2<=-4823/11250<-2/5 at xi=1/6",
            "samplers": "every fixed real nonzero g in C_c^infinity(-1,1); the same family works for all such g",
            "divergence": "integral g^2*T_ll <-(3/10)*hbar*R^2*gamma0^2*||g||^2 tends to minus infinity",
            "same_one_sided_cap": "w_R<=Phi_*^2 for every fixed Phi_*^2>=0 on the entire sampled null segment",
            "finite_bound_depending_only_on_sampler_and_nonnegative_Wick_upper_cap_exists": False,
            "state_cap_holds_everywhere_or_is_two_sided": False,
            "arbitrary_fixed_affine_length": "dilated normalized mode tau*h_R(tau*k), field and stress weights tau^-1,tau^-4",
        },
        "energy_and_ANEC": {
            "actual_total_energy": "hbar*(R^2*Jp*I0+J2*I2)/(6*Iq*J2)>0 and finite for each finite R",
            "improvement_spatial_boundary_integrates_to_zero": True,
            "family_has_common_finite_energy_or_physical_momentum_budget": False,
            "positive_frequency_support_kills_anomalous_complete_line_term": True,
            "actual_complete_line_ANEC": "integral T_ll ds=(4*pi/3)*A_R^2*I2>0 at every finite R",
            "new_general_ANEC_theorem_or_negative_ANEC_claimed": False,
            "limit_state_or_interchange_of_divergent_limits_claimed": False,
        },
        "verification_boundary": "Native exact mass-shell/Jacobian and normalized continuum-mode identities, actual positive pure oscillator covariance, parent physical stress contraction including improvement, continuous Taylor/mollifier jet bounds, independent Fraction constants and actual finite-energy/positive-ANEC formulas. Hadamard, Fock implementation, Fourier support and whole-segment proofs are written and source-pinned, not proof-assistant formalized or independently peer reviewed.",
        "not_established": [
            "Two-sided Wick, total-energy or physical-momentum capped null bounds",
            "Transverse/double smearing, entropy or alternative null focusing inequalities",
            "A self-consistent SEE counterexample, actual gravitational incompleteness or bounce UV verdict",
            "General ANEC beyond this explicit family, arbitrary realistic/interacting fields or original P8 closure",
        ],
        "verdict": "An actual normalized smooth one-mode squeezed Hadamard family has unboundedly negative compact null energy while its relative Wick square is nonpositive on the whole sampled segment, including minimal and conformal couplings. A one-sided nonnegative field-strength upper cap cannot rescue a finite null-line lower bound. Each state has finite energy and strictly positive complete-line ANEC; no common UV/resource or two-sided cap is maintained. This rules out the direct null extension of the timelike state-cap strategy, not all null singularity routes or original P8.",
    }
    assert tuple(result) == REPORT_KEYS
    return result


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete null state-cap certificate differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), result)
        print(
            "P8 A.21 actual null one-sided-state-cap obstruction replay passed; general null routes and original P8 OPEN"
        )
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
