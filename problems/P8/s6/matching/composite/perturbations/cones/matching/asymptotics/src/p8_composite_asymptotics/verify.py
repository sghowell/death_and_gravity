"""Read-only scaled-limit identities and fixed-interval theorem replay."""

import argparse
import hashlib
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

import sympy as sp
from p8_composite_light import verify as prior

from . import independent, limiting, scaled

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"composite-asymptotics.json"
PRIOR_SHA = "36c899d2b82957f2fc83df8e7bc66582cbf77829fbed8b64724d2d39c2086e48"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(value):
    if isinstance(value, dict):
        return {str(key): encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    return str(value)


@cache
def prior_checks():
    if sha(prior.REPORT) != PRIOR_SHA:
        raise ValueError("Pinned S6.9 light-mode report changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return PRIOR_SHA


def controls():
    data = limiting.general_jets()
    at = {data["v"]: 2, data["a"]: 11}
    true_acceleration = data["heavy_growth_coefficient"].subs(at)
    omitted_acceleration = -data["mass_squared"].subs(at)
    backward_scale = 7*sp.Rational(9, 10)+5/sp.Rational(9, 10)-11
    k = sp.Symbol("kappa", positive=True)
    e = sp.Symbol("e", nonnegative=True)
    source = (e*k-1)/((1+e)*sp.sqrt(k))
    return {
        "true_initial_heavy_acceleration": true_acceleration,
        "omitted_normalization_false_acceleration": omitted_acceleration,
        "lost_normalization_acceleration": true_acceleration-omitted_acceleration,
        "backward_z_below_one_margin": 1-backward_scale,
        "nonzero_normalized_source_limit": source.subs(e, 0),
        "small_metric_amplitude_absolute_projection_limit": sp.limit(e*source, e, 0),
        "zero_heavy_data_unforced_solution": sp.S.Zero,
    }


@cache
def build_report():
    previous = prior_checks()
    residuals = {
        "full_scaled_positive_root_CD_ODE": scaled.checks(),
        "exact_forward_limiting_flow_and_comparison": limiting.checks(),
        "complete_canonical_coefficient_limits": limiting.canonical_limit_checks(),
        "physical_projection_time_and_tidal_jets": limiting.physical_projection_checks(),
    }
    if any(sp.simplify(value) != 0 for group in residuals.values() for value in group.values()):
        raise ValueError("A scaled-flow or actual limiting response identity failed")
    rational = independent.checks()
    omissions = controls()
    if not (
        omissions["true_initial_heavy_acceleration"] == sp.Rational(55, 2)
        and omissions["omitted_normalization_false_acceleration"] == -2
        and omissions["lost_normalization_acceleration"] == sp.Rational(59, 2)
        and omissions["backward_z_below_one_margin"] == sp.Rational(13, 90)
        and omissions["small_metric_amplitude_absolute_projection_limit"] == 0
        and Fraction(140, 19) < 9
    ):
        raise ValueError("A direction, amplitude or normalization control failed")
    general = limiting.general_jets()
    trajectory = limiting.solution()
    sources = sorted(ROOT.glob("src/p8_composite_asymptotics/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1,
        "claim": "P8-S6.10.COMPOSITE",
        "date": "2026-09-06",
        "status": "EXACT_LARGE_ASYMMETRY_LIMITING_HOMOGENEOUS_RESPONSE_AND_FIXED_INTERVAL_CONVERGENCE; GENERAL_MATCHING_AND_P8_OPEN",
        "prior_S6_9_composite_light_sha256": previous,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
        "formulation": "FORMULATION.md",
        "written_proof": "notes/proof.md",
        "source_audit": "notes/sources.md",
        "exact_residuals": {name: dict.fromkeys(values, "0") for name, values in residuals.items()},
        "independent_Fraction_replay": rational,
        "omission_amplitude_and_direction_controls": encode(omissions),
        "domain": {
            "action": "Constant G=F=M^2,m4=M^2*m^2,alpha=beta=1,beta_n=(0,0,1,0,0), canonical shared composite matter",
            "physical_frame": "Pinned +--- composite g_eff, physical u=m*T, Ae(0)=1; not original CD/M1 frame",
            "finite_parameter_family": "epsilon=1/Y>0,y0=Y,rhobar0=9/Y^2,h=11epsilon*u/(1+11epsilon*u^2/4)",
            "CD_duration": "m*tau=2/sqrt(11epsilon)",
            "scaled_variables": "e=1/y=epsilon*z,eta=rhobar/e^2,j=h/e,a_current=hprime/e",
            "finite_background_existence": "Sufficiently small positive epsilon depending on each fixed forward interval [0,L]",
            "matter": "Positive canonical null source with a locally reconstructed potential for each finite member; density/potential may be negative",
            "mode": "Exact fixed-charge k=0 relative tensor, J=0,l0=0,H0=1,Hprime0=0; unit normalized linear column",
            "time_direction": "u>=0 only for the explicit orbit bounds; not a two-sided bounce-duration theorem",
        },
        "exact_limiting_trajectory": encode(trajectory),
        "general_leading_jets": encode(general),
        "forward_interval_bounds": {
            "z": "z>=1, zprime>0",
            "R": "(zprime-11u)prime=z>=1 and initial numerator=2 imply R>0",
            "v_squared": "v^2=1+22/z-19/z^2, 1<=v^2<=140/19<9",
            "upper_square": "140/19-v^2=(11z-19)^2/(19z^2)",
            "heavy_coefficient": "D_H=N_R-m_alg^2=3(a+1)^2/(4v^2)+(-v^2+4v-2)/4",
            "positive_margin": "D_H-1/4=3(a+1)^2/(4v^2)+(v-1)*(3-v)/4>=0",
        },
        "actual_limiting_homogeneous_response": {
            "equation": "Hsecond=D_H(u)*H, H0=1,Hprime0=0; l=0 in the exact limit",
            "comparison": "H(u)>=cosh(u/2) for every u>=0",
            "positive_retarded_identity": "H-cosh(u/2)=integral_0^u 2sinh((u-s)/2)*(D_H(s)-1/4)*H(s) ds",
            "state_control": "Zero heavy data give the exact zero homogeneous solution; growth is not compulsory for all states",
            "scope": "Full time-dependent limiting k=0 equation, not a frozen diagonal or a finite-k spectral statement",
        },
        "physical_projection": {
            "exact": "gamma_eff=(l+Rsrc*H)/f_sum, Rsrc=(e*kappa-1)/[(1+e)sqrt(kappa)]",
            "limit": "f_sum=M/2,Rsrc=-sqrt(v),gamma_eff=-(2/M)sqrt(v)*H",
            "fixed_gauge_amplitude_bound": "abs(gamma_eff)>=2cosh(u/2)/M",
            "constant_gauge_invariant_difference": "abs(gamma_eff(u)-gamma_eff(0))>=2max{0,cosh(u/2)-sqrt(2)}/M",
            "initial_time_jet": "gamma_eff_prime(0)=-4sqrt(2)/M in u time",
            "initial_tidal_jet": "gamma_eff_second(0)=-20sqrt(2)/M in u time, nonzero physical linearized tidal curvature",
        },
        "fixed_interval_convergence": {
            "quantifiers": "For each fixed finite L and tolerance, sufficiently small positive epsilon gives regular finite backgrounds and uniform convergence on [0,L]",
            "background": "Analytic scaled vector field on a compact positive-root tube; error <=B_L*epsilon*(exp(K_L*L)-1)/K_L",
            "canonical": "Regularized kinetic/coefficient jets are analytic on the same tube despite f_relative=O(epsilon)",
            "mode": "Matched J=0 canonical data give (H,Hprime,l)->(H_limit,Hprime_limit,0), O_L(epsilon)",
            "physical": "Composite projection and required time derivatives converge on the same fixed interval",
            "not_computed": "Numerical tube/propagator bounds, epsilon_*(L), or an estimate on growing intervals",
        },
        "singular_limit_and_amplitude_guard": {
            "individual_geometry": "N_g and a_g vanish as epsilon->0; endpoint is not a regular two-metric parent solution",
            "scalar_clock": "Positive at finite epsilon but degenerates in the limit; no uniform reconstructed-potential derivative bounds",
            "normalized_column": "H0=1 implies relative metric delta=H/f_relative=O(1/epsilon), not uniform small metric amplitude",
            "small_metric_family": "Multiply the linear column by epsilon*epsilon_lin on a fixed interval; relative metrics remain small and absolute composite amplitude tends to0",
            "surviving_statement": "Normalized linear growth ratios and tidal derivatives, not finite nonlinear instability or unsuppressed absolute response at fixed small metric amplitude",
        },
        "general_fixed_scaling_screen": "At a centre with fixed eta0>0,a0>=0, v0=sqrt(1+eta0/3)>1 and N_R-(3/4)m_alg^2=[3(a0+1)^2+v0^2(3v0-2)]/(4v0^2)>0",
        "verification_boundary": "Written positivity/comparison/compact-Gronwall proofs plus exact source/ODE/canonical/Fraction identities and independent audits; not Lean formalized, no numerical finite-Y threshold",
        "not_established": [
            "A regular two-metric geometry at epsilon=0",
            "Finite-Y all-time instability, a nonlinear instability, or growth for every initial state",
            "A uniform finite-amplitude physical two-metric perturbation family for H0=1",
            "A bound on a CD-duration or logarithmically growing window, or an interchange of Y and time limits",
            "A finite-k gap, dispersion relation, vector/scalar health, cutoff, ghost or ultraviolet verdict",
            "A general exclusion of nonadiabatic finite-band light-only reductions or other HR states/actions",
            "A globally free old M1 matter action or original CD/M1 physical-frame matching",
            "Completion of S6, P8(b), or P8",
        ],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Composite asymptotic-response certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.10.COMPOSITE: exact limiting response and fixed-interval replay passed; general matching and P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
