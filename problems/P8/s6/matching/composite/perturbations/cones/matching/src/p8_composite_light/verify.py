"""Read-only exact rolling two-TT algebra and specified hierarchy replay."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_composite_cones import verify as prior

from . import background, bounds, canonical, independent

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"composite-light.json"
PRIOR_SHA = "89145c7a91981f72941a30e902b7e83aa231f4f4c97440c8e6f454ba339c340a"


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
        raise ValueError("Pinned S6.8 composite-cone report changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return PRIOR_SHA


def domain_controls():
    result = {}
    for label, value in (("ordinary_float", 0.1), ("symbolic_float", sp.Float("0.1")),
                         ("infinite", sp.oo), ("not_a_number", sp.nan)):
        try:
            bounds.linear_root_enclosure(value)
        except (TypeError, ValueError, sp.PolynomialError):
            result[label] = "rejected"
        else:
            raise ValueError("An inexact or nonfinite enclosure was accepted")
    t = sp.Symbol("T", real=True)
    try:
        canonical.mass_led_pair(t, t, 1, 0, 0, 0)
    except ValueError:
        result["zero_mass_inverse_representative"] = "rejected"
    else:
        raise ValueError("The inverse-mass representative accepted zero mass")
    if canonical.generic_action()["mass_squared"].is_positive is True:
        raise ValueError("The exact canonical identity silently imposed positive mass")
    result["general_mass_domain"] = "real; zero and negative allowed in exact identities"
    return result


@cache
def build_report():
    previous = prior_checks()
    residuals = {
        "physical_clock_and_weighted_fields": canonical.physical_checks(),
        "full_action_boundary_Eulers_and_Hamiltonian": canonical.generic_action_checks(),
        "retarded_and_inverse_mass_equation_identities": canonical.operator_checks(),
        "exact_fixed_charge_zero_k_reduction": canonical.routh_checks(),
        "actual_full_reconstruction_ODE_jets": background.checks(),
    }
    if any(sp.simplify(value) != 0 for group in residuals.values() for value in group.values()):
        raise ValueError("An exact rolling two-TT identity failed")
    rational = independent.checks()
    enclosures = bounds.build()
    omissions = canonical.controls()
    if omissions != {
        "omitted_rotating_weight_heavy_forcing": sp.S(2),
        "omitted_homogeneous_heavy_data_light_forcing": sp.S(-2),
        "omitted_nonconstant_mass_derivatives_residual": sp.S(-12),
        "omitted_time_boundary_changes_momentum": sp.S(-1),
    }:
        raise ValueError("A time-dependent or state omission control failed")
    physical = canonical.physical_coefficients()
    centres = {str(y): background.initial_jets(y) for y in (1, 2)}
    at = {background.A: sp.Rational(4, 25)}
    sample_keys = ("mass_squared", "c_light_squared", "omega", "theta_sum",
                   "N_sum", "N_relative", "Routh_frequency_squared")
    samples = {y: {key: sp.simplify(data[key].subs(at)) for key in sample_keys}
               for y, data in centres.items()}
    if samples["1"]["Routh_frequency_squared"] != sp.Rational(153, 100):
        raise ValueError("The positive zero-mass Routh countercontrol failed")
    sources = sorted(ROOT.glob("src/p8_composite_light/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1,
        "claim": "P8-S6.9.COMPOSITE",
        "date": "2026-09-06",
        "status": "EXACT_ROLLING_TWO_TT_AND_SPECIFIED_MASS_LED_HIERARCHY_SCREENS; GENERAL_LIGHT_ONLY_MATCHING_AND_P8_OPEN",
        "prior_S6_8_composite_cones_sha256": previous,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
        "formulation": "FORMULATION.md",
        "written_proof": "notes/proof.md",
        "source_audit": "notes/sources.md",
        "exact_residuals": {name: dict.fromkeys(values, "0") for name, values in residuals.items()},
        "independent_Fraction_replay": rational,
        "omission_controls": encode(omissions),
        "domain_controls": domain_controls(),
        "domain": {
            "action": "Constant G,F,m4>0, constant real beta_n, constant alpha,beta>0; canonical matter shared through the composite metric",
            "geometry": "Sufficiently smooth regular common-flat positive-root geometry, positive finite lapses/scales; no denominator degeneracy",
            "physical_frame": "+--- g_eff=g(alpha I+beta sqrt(g^-1 f))^2, Ae=a*r, Ne=Ng*s, dT=Ne*dt",
            "TT_scope": "Two independent TT polarizations, identical two-field action per real contraction; no scalar/vector elimination",
            "wave_number": "q=k/Ae for fixed comoving k; exact fixed-charge reduction is separately k=0 only",
            "general_mass": "Real m_alg^2, including zero or negative; nonzero interval needed only for the formal inverse-mass representative",
            "screen_model": "G=F=M^2,m4=M^2*m^2,alpha=beta=1,beta_n=(0,0,1,0,0)",
            "screen_data": "rho/(M^2*m^2)=1/2,y0=1 or 2,h0=hsecond0=0,hprime0=A>=0",
            "CD_duration": "u=m*T,h=He/m,A=4/(m*tau)^2 for finite CD duration; A=0 is included only as an algebraic endpoint",
            "matter_matching": "Local analytically reconstructed canonical potential; not one globally free old M1 action across the parameter family",
        },
        "physical_coefficients": encode({key: physical[key] for key in
                                         ("K1", "K2", "K_relative", "relative_stiffness_T",
                                          "c_light_squared", "c_heavy_squared", "D", "mass_squared")}),
        "canonical_map_and_boundary": {
            "map": "l=f_sum*(w1*h+w2*gamma), H=f_relative*(gamma-h), f_sum=sqrt(2K_sum), f_relative=sqrt(2K_relative)",
            "inverse": "h=l/f_sum-w2*H/f_relative; gamma=l/f_sum+w1*H/f_relative",
            "rotation": "omega=sqrt(w1*w2)*dT log sqrt(K2/K1)=f_sum*w2prime/(2*f_relative)",
            "boundary": "Lpre=Lpost-dT[(theta_sum*l^2+theta_relative*H^2)/2]",
            "post_momenta": "p_l=lprime-2omega*H, p_H=Hprime",
            "Euler_diagonals": "V_LL=q^2*cL^2-N_sum; V_HH=m_alg^2+q^2*cH^2-N_relative-4omega^2",
            "adjoint": "B=2omega*(dT-theta_sum)+q^2D; Bstar=-2omega*dT-2omegaprime-2omega*theta_sum+q^2D",
        },
        "retarded_and_local_representative": {
            "retarded_solution": "H=H_hom-G_R*B*l for specified heavy initial data and zero-initial-response retarded G_R",
            "exact_light_equation": "(D_L-Bstar*G_R*B)l+Bstar*H_hom=0",
            "causal_action_guard": "The retarded inverse is not substituted as a symmetric single-copy effective action",
            "local_pair": "H0=-B*l/m_alg^2, D_L*l-Bstar*(B*l/m_alg^2)=0, only on nonzero mass interval",
            "exact_heavy_residual": "-(dT^2+V_HH-m_alg^2)[B*l/m_alg^2]",
            "formal_local_action_correction": "+(B*l)^2/(2*m_alg^2), with complete coefficient derivatives",
            "conditional_approximation": "A specified full propagator C_I gives phase error <=C_I*(initial phase error+integral abs(R_H)); no advantageous C_I or small bound supplied here",
            "required_research_inputs": [
                "Uniform positive gap and derivative estimates for a proposed mass-led expansion",
                "Chosen finite band, light derivative class and bounded mixing/background rates",
                "Actual heavy/full Green bound and matched or suppressed homogeneous heavy data",
                "Physical inverse-map bounds and omitted scalar/vector/cutoff/matching analysis",
            ],
        },
        "fixed_charge_zero_momentum": {
            "charge": "J=2K_sum*(uprime-w2prime*delta)=f_sum*(p_l-theta_sum*l), Jprime=0",
            "exact_equation": "Hsecond+(m_alg^2-N_relative)*H=-w2prime*J/f_relative",
            "countercontrol": "At y0=1,A=4/25: m_alg^2=0 but Omega_Routh^2/m^2=153/100>0",
            "scope": "Fixed-charge k=0 equation only; neither full finite-k spectrum nor a positive retarded Green bound",
        },
        "actual_centre_jets": encode(centres),
        "samples_at_m_tau_5": encode(samples),
        "strict_hierarchy_enclosures": encode(enclosures),
        "units": {
            "physical_time": "T is composite proper time; u=m*T only in centre-jet evaluation",
            "first_order": "theta,omega,He have dimension mass; centre jets report them divided by m",
            "second_order": "N,m_alg^2,Omega_Routh^2,Hdot_eff have dimension mass^2; centre jets report them divided by m^2",
            "comparison": "Hdot_eff/m^2=A; all hierarchy ratios use the same physical clock",
        },
        "proved_screen": [
            "y0=1: zero algebraic mass prevents the mass-led inverse-mass representative through the centre",
            "y0=2: 0<m_alg^2/m^2<9/200 and 0<c_locked^2<1",
            "y0=2,all A>=0: max{Hdot_eff,omega^2}/m_alg^2>50/9",
            "Both specified states: N_sum>0 and Omega_Routh^2<3N_sum for A>=0",
            "Changing m*tau cannot produce the stated parametric mass-led/slow-background hierarchy at these specified data",
        ],
        "verification_boundary": "Written exact action/ODE/retarded-identity and rational inequality proof plus separate Fraction and independent covariant audits; not a computer-formalized analytic propagator theorem",
        "not_established": [
            "A finite-k physical spectral gap, eigenmode evolution, growth rate, or positive heavy Green bound",
            "The absence of all canonical frequencies from zero mu or a negative diagonal canonical potential",
            "Failure of every nonadiabatic finite-window light-only reduction or every other canonical split",
            "A statement about all initial ratios, densities, HR coefficients, or parent states",
            "A general composite/bimetric/ultraviolet exclusion or application of the all-mode S6.8 contract to a light-only theory",
            "A scalar/vector health result, cutoff, quantum cone, state, positivity amplitude, or global vacuum",
            "Matching to the original frozen CD/M1 physical metric and free matter action",
            "Completion of S6, P8(b), or P8",
        ],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Composite light-mode certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.9.COMPOSITE: exact rolling TT and specified hierarchy replay passed; general light-only matching and P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
