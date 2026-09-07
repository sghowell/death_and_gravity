"""Read-only exact finite-positive-parameter source-response certificate replay."""

import argparse
import hashlib
import json
from fractions import Fraction as Q
from functools import cache
from pathlib import Path

import sympy as sp
from p8_composite_asymptotics import verify as prior

from . import bounds, independent, model, response

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"composite-response.json"
PRIOR_SHA = "bfb0b74a2d50a508a5e5fc346b59a31958c9f9a8ccfad621d4caa34883e90eef"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@cache
def prior_checks():
    if sha(prior.REPORT) != PRIOR_SHA:
        raise ValueError("Pinned S6.10 composite-asymptotic report changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return PRIOR_SHA


def controls():
    eps = sp.Symbol("epsilon", nonnegative=True)
    data = {
        "finite_family_initial_R_parameter_derivative": sp.diff(model.initial_root(eps), eps).subs(eps, 0),
        "wrong_constant_initial_R_changes_eta_derivative": sp.diff(9/(1+eps)**3, eps).subs(eps, 0),
        "limiting_retarded_physical_kernel_initial_derivative": Q(3),
        "locked_kernel_initial_derivative": Q(1),
        "omitted_relative_source_output_contribution": Q(2),
        "larger_band_cannot_use_positive_kernel_bound": Q(1, 4)-Q(140, 19)/4,
        "opposite_source_pulse_linear_moment": -Q(1, 60),
        "missing_locked_kinetic_at_epsilon_1_over_10000_R1": Q(1, 10001**2),
        "zero_source_zero_data_response": Q(0),
    }
    expected = (sp.Rational(9, 4), -27, Q(3), Q(1), Q(2), -Q(121, 76),
                -Q(1, 60), Q(1, 10001**2), Q(0))
    if tuple(data.values()) != expected:
        raise ValueError("A source, band, initial-family or locked-action control failed")
    return {key: str(value) for key, value in data.items()}


@cache
def build_report():
    previous = prior_checks()
    residuals = {"exact_background_flux_and_locked_action": model.checks(),
                 "limiting_normalization_and_pulse": response.checks(),
                 "physical_clock_probe_and_metric_source": response.source_checks()}
    if any(sp.simplify(value) != 0 for group in residuals.values() for value in group.values()):
        raise ValueError("An exact background, source or limiting response identity failed")
    interval = bounds.build()
    comparisons = response.rational_checks()
    rational = independent.checks()
    omissions = controls()
    sources = sorted(ROOT.glob("src/p8_composite_response/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1,
        "claim": "P8-S6.11.COMPOSITE",
        "date": "2026-09-06",
        "status": "EXPLICIT_FINITE_PARAMETER_PHYSICAL_TT_SOURCE_RESPONSE_SCREEN_FOR_LOCKED_ACTION; GENERAL_MATCHING_AND_P8_OPEN",
        "prior_S6_10_composite_asymptotics_sha256": previous,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
        "formulation": "FORMULATION.md",
        "written_proof": "notes/proof.md",
        "source_audit": "notes/sources.md",
        "exact_residuals": {name: dict.fromkeys(group, "0") for name, group in residuals.items()},
        "independent_Fraction_replay": rational,
        "source_band_initial_data_and_omission_controls": omissions,
        "outward_rational_finite_parameter_bounds": interval,
        "exact_comparison_and_tube_margins": comparisons,
        "domain": {
            "constant_model": "G=F=M²>0,m4=M²m²,m>0,alpha=beta=1,beta_n=(0,0,1,0,0)",
            "frame": "Positive-root +--- composite effective metric, not original CD/M1 matter frame",
            "internal_matter": "Shared canonical scalar with the pinned family's reconstructed potential; positive null source, density/potential need not be positive",
            "family": "epsilon=1/Y>0,z0=1,eta0=9,R0=sqrt(1+3(1+epsilon)^3),Ae=(1+11epsilon*u²/4)²",
            "time": "u=mT in [0,1], physical composite proper time T",
            "spatial_band": "kbar=kcom/m in [1/16,1/8]; physical effective momentum m*kbar/Ae(u)",
            "parameter": interval["epsilon_range"],
            "small_metric_amplitude": interval["sigma_amplitude_range"],
            "retarded_data": "Both TT metric fields and their fluxes vanish initially, and the locked candidate has the same zero data",
        },
        "physical_probe": {
            "polarization": "One real TT polarization eij eij=1, eii=0, ki eij=0; both polarizations obey identical equations",
            "stress_convention": "delta Sprobe=(1/2) integral sqrt|g_eff| T^ij delta g_eff,ij; T^ij=Pi eij/Ae²,T^0mu=0,delta g_eff,ij=-Ae²Gamma eij",
            "background_conservation": "Time divergence is proportional to spatial trace and spatial divergence to ki eij; both vanish for arbitrary C1 Pi(T)",
            "pulse": "p(u)=u²(1-u)² on [0,1], extended by zero, C1",
            "literal_u_action_source": "J_u=m M² epsilon sigma p/4",
            "action_divided_by_m_source": "Jbar=J_u/m=M² epsilon sigma p/4",
            "physical_stress": "Pi=-2m J_u/Ae³=-M²m² epsilon sigma p/(2Ae³)",
            "scope": "External infinitesimal conserved probe, not supplied by the internal scalar background; no microscopic probe or temporal-bandlimit claim",
        },
        "exact_regular_flux_problem": {
            "variables": "g=h_g/sigma,f=h_f/(epsilon sigma),Pg=Ag*gprime,Pf=Af*fprime",
            "output": "R=Gamma/(epsilon sigma)=(zg+f)/(1+e)",
            "locked": "g=epsilon L,f=L; (C_lock Lprime)prime+V_lock L=p, C_lock=Af+epsilon²Ag,V_lock=Vf+epsilon²Vg",
            "regularity": "All required kinetic/constraint denominators positive on certified finite member; no division by relative mass or Hubble rate",
            "constraints": "The external TT stress carries no linear scalar/vector source; full background constraints replay the pinned reconstruction",
        },
        "limiting_causal_response": {
            "light": "fsecond+q*f=p",
            "relative": "chi=(z/sqrt(v))*g; chisecond+(v²q-D_H)*chi=sqrt(v)*p",
            "physical_output": "R0=f+sqrt(v)*chi",
            "kernel": "K0(u,s)=G_L(u,s)+sqrt(v(u))*G_H(u,s)*sqrt(v(s))",
            "comparison": "D_H-v²q>=41/304>1/8; G_H>=sqrt(8)sinh((u-s)/sqrt(8))>=u-s; (383/384)(u-s)<=G_L<=u-s",
            "pulse_at_u1": "I=1/60, (383/384)I<=L0<=I, relative contribution H0>=I, R0=L0+H0",
            "causality": "Explicit zero-data initial value problem, no replacement of retarded inverse by a symmetric single-copy action; other homogeneous data would add a response",
        },
        "finite_parameter_theorem": {
            "uniform_output_errors": "Both |R_epsilon-R0| and |L_epsilon-L0| <=2^-10<1/600 on the whole fixed interval and band",
            "u1_absolute_mismatch": "R_epsilon-L_epsilon>=1/75",
            "u1_positive_outputs": "Both R_epsilon and L_epsilon strictly positive for the positive pulse",
            "u1_relative_mismatch": "(R_epsilon-L_epsilon)/R_epsilon>=1/3; (2/3)R_epsilon-L_epsilon>=1/360",
            "small_metric": "Each individual metric polarization amplitude <=1/128<1/100 under the sigma gate",
            "absolute_vs_normalized": "Physical probe and composite response vanish as epsilon*sigma; the source-normalized relative error does not",
            "constants": "Explicit sufficient unoptimized exact rational gates, not floating-point propagation or a physical cutoff estimate",
        },
        "prepared_source_free_positive_control": {
            "data": "p=0; h_g=h_f=epsilon*sigma*l0 and equal u derivatives epsilon*sigma*v0, |l0|,|v0|<=1; exact H0=Hprime0=0",
            "full_initial_flux": "(epsilon*l0,epsilon*Ag*v0,l0,Af*v0,l0,C_lock*v0)",
            "limit": "g=0,f=L and both physical outputs agree",
            "finite_error": "|R_epsilon-L_epsilon|<=2^187104*epsilon<1/600 under the same parameter gate",
            "scope": "Prepared source-free light data, not universal matching under arbitrary physical sources or relative initial data",
        },
        "verification_boundary": "Written continuation/Volterra/Gronwall proof, exact symbolic and separate Fraction algebra, outward interval AD and source-hashed immutable replay; not proof-assistant formalized",
        "not_established": [
            "Failure of every second-order EFT or every nonadiabatic light-only reduction without its own source, memory, locality and error contract",
            "A regular two-metric geometry at epsilon zero or a uniform finite-nonlinear-amplitude theorem",
            "Uniform small parent proper-time curvature, an EFT cutoff hierarchy, or a temporal-frequency-limited probe",
            "A finite-k spectral, scalar/vector-health, ghost, quantum-state or ultraviolet verdict",
            "Uniformity on a growing logY interval or a full CD-duration window",
            "Original-frame CD/M1 matching, globally free internal M1 matter, or completion of S6/P8",
        ],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Composite finite-source-response certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.11.COMPOSITE: explicit finite-parameter physical source-response replay passed; general matching and P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
