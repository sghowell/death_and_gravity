"""Read-only source-hashed exact actual-auxiliary-metric cone replay."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_trimetric import verify as prior

from . import background, bounds, independent, rolling, tensor

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"actual-tensor-cones.json"
P8 = next(path for path in ROOT.parents if path.name == "P8")
PRIOR_SHA = "062ea71b4fcb139af3ea727f013eeeb3a340b6e37207cb872e372eacce09631e"
P8_FORMULATION_SHA = "73dd18f7a4b56a0205f9fbc6a4b09b213a3c280b3f0274414a255635b04e2b09"
CD_SHA = "caf8c8e688a7565b9d00f921c099a28da00f97522ed26ad182a8227eb80cd4dd"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@cache
def prior_checks():
    pins = {"S6_13_actual_auxiliary_parent": (prior.REPORT, PRIOR_SHA),
            "original_B_curvature_and_physical_cone_contract": (P8/"FORMULATION.md", P8_FORMULATION_SHA),
            "conditional_original_CD_M1_free_chi": (P8/"certificates"/"witness-CD_matter.json", CD_SHA)}
    if any(sha(path) != expected for path, expected in pins.values()):
        raise ValueError("A pinned action, convention or original M1 input changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {key: expected for key, (_, expected) in pins.items()}


def controls():
    b, t, c = background.controls(), tensor.controls(), bounds.controls()
    if tuple(b.values()) != (0, 0, sp.Rational(1, 2), 1, 2, 0, 0, 0, 4, 0):
        raise ValueError("A mixed-sign auxiliary/full-vacuum or nonrolling control failed")
    td = tensor.derive()
    if tuple(t.values()) != (48, 1, 0, (td["G"]-td["F"])/4, 0):
        raise ValueError("A physical clock, source contact or symmetry control failed")
    if tuple(c.values()) != (6, 48, sp.Rational(24, 7), sp.Rational(3, 56), sp.Rational(95, 28), 0):
        raise ValueError("A physical coefficient-error threshold control failed")
    return {"background_hypotheses": {key: str(value) for key, value in b.items()},
            "tensor_source_clock_and_symmetry": {key: str(value) for key, value in t.items()},
            "conditional_error_budgets": {key: str(value) for key, value in c.items()}}


def clock_bridge():
    """Compare genuinely separate Fraction measure/derivative pullbacks."""
    d = tensor.derive()
    bridge = []
    for fixture in independent.clock_fixtures():
        point = {d[key]: sp.Rational(fixture[key].numerator, fixture[key].denominator)
                 for key in ("A", "N", "G")}
        residuals = {}
        for key in ("G_T_h", "F_T_h", "common_speed_squared"):
            oracle = sp.Rational(fixture[key].numerator, fixture[key].denominator)
            if sp.cancel(d[key].subs(point)-oracle) != 0:
                raise ValueError("Primary versus independent physical-clock pullback failed")
            residuals[key] = "0"
        bridge.append({"inputs": {key: str(fixture[key]) for key in ("A", "N", "G")},
                       "primary_minus_independent": residuals})
    return bridge


@cache
def build_report():
    previous = prior_checks()
    residuals = {"full_covariant_background_cones": background.checks(),
                 "literal_TT_source_and_physical_clock": tensor.checks(),
                 "actual_free_scalar_local_rolling_solution": rolling.checks(),
                 "conditional_coefficient_error_thresholds": bounds.checks()}
    if any(sp.simplify(value) != 0 for group in residuals.values() for value in group.values()):
        raise ValueError("A background, TT, actual solution or bound identity failed")
    td, rd = tensor.derive(), rolling.derive()
    if any(td[key].is_positive is not True for key in ("G_T_h", "F_T_h", "relative_algebraic_mass_h_squared")):
        raise ValueError("The exact tensor coefficient domains failed")
    if any(rd[key].is_positive is not True for key in ("A", "N", "rho", "H_r_squared", "D_h_psi")):
        raise ValueError("The actual rolling solution positivity domain failed")
    rational, omissions = independent.checks(), controls()
    sources = sorted(ROOT.glob("src/p8_trimetric_cones/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1, "claim": "P8-S6.14.AUXILIARY", "date": "2026-09-06",
        "status": "POSITIVE_LINK_FULL_TENSOR_CONE_OBSTRUCTION; SYMMETRIC_PHYSICAL_COMMON_CHANNEL_EXCEEDS_MATTER_CONE; P8_OPEN",
        "prior_context_sha256": previous,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
        "formulation": "FORMULATION.md", "written_proof": "notes/proof.md", "source_audit": "notes/sources.md",
        "exact_residuals": {name: dict.fromkeys(group, "0") for name, group in residuals.items()},
        "independent_Fraction_replay": rational,
        "independent_physical_clock_bridge": clock_bridge(), "checked_controls": omissions,
        "action_and_domain": {
            "action": "The pinned S6.13 auxiliary action, now allowing explicit -2beta4g det(e)-2beta4f det(v); no u EH or other derivative terms",
            "constants": "G,F,p_g,p_f,epsilon>0; B,beta4g,beta4f constant real; positivity of links is an explicit premise, not inferred from vacuum mass",
            "actual_physical_metric": "h=u^T eta u, not its vacuum approximation and not either Einstein metric",
            "B_curvature": "+---; R_B=-6(Hdot+2H²), Einstein_B00=+3H², EH=-G R_B/2; not the opposite A/FK curvature convention",
            "geometry": "Actual regular common spatially flat smooth solution, all six lapse/scale entries positive, parent symmetrization and root domain",
            "matter_for_background_identity": "Homogeneous isotropic density rho and pressure p; n_h=epsilon(rho+p)",
            "matter_for_tensor_action": "Homogeneous positive-field-metric canonical scalar(s) on h, with no independent TT anisotropic response; isotropy alone is insufficient for a general medium",
            "physical_normalization": "epsilon multiplies bare matter; fixed canonically normalized free-chi stress is not suppressed by changing epsilon",
        },
        "general_full_parent_cone_theorem": {
            "weights_and_speeds": "P_g=p_g a_e/a_u,P_f=p_f a_v/a_u,c_e=a_u n_e/(n_u a_e),c_f=a_u n_v/(n_u a_v)",
            "full_u_equations": "B+3(P_g+P_f)=-epsilon rho/2; B+P_g(c_e+2)+P_f(c_f+2)=epsilon p/2",
            "undivided_identity": "P_g(c_e-1)+P_f(c_f-1)=n_h/2, with no Hubble, B or source-expansion division",
            "positive_weight_bound": "max(c_e,c_f)>=1+n_h/[2(P_g+P_f)]>1 if n_h>0",
            "scope": "At least one full EH tensor principal cone, not both individually and not an asymmetric finite-band light-only matching exclusion",
            "quadratic_light_countercontrol": "K=diag(100,1),gradient=diag(25,4),positive heavy relative spring: high squared cones 1/4,4 but locked light squared speed 29/101; not claimed a parent background",
        },
        "exact_auxiliary_TT_and_external_probe": {
            "polarization": "tr(E)=0,tr(E²)=1, spatial coframes exp(gamma_i E/2)",
            "algebraic_density": "L_alg/sqrt|h|=-[P_g(gamma_e-gamma_u)²+P_f(gamma_v-gamma_u)²]/4+j gamma_u",
            "source_normalization": "S_j=int dT d³x a_h³ j gamma_u; for canonical variational stress deltaS=-1/2 int sqrt|h| T^munu deltah_munu, j=Pi_TT/2 (opposite plus-variation stress uses j=-Pi_TT/2)",
            "source_scope": "External infinitesimal physical TT probe, separately background-conserved for any smooth time profile, not the internal canonical scalar and no probe positivity assumption",
            "auxiliary_map": "gamma_u=(P_g gamma_e+P_f gamma_v+2j)/(P_g+P_f)",
            "stationary_action": "-P_gP_f(gamma_e-gamma_v)²/[4(P_g+P_f)]+j(P_g gamma_e+P_f gamma_v)/(P_g+P_f)+j²/(P_g+P_f)",
            "boundary": "The local source contact is retained and does not cancel a propagating common response; no scalar auxiliary inverse is inferred from this nonzero TT denominator",
        },
        "symmetric_physical_common_channel": {
            "additional_hypotheses": "G=F,p_g=p_f=q>0,beta4g=beta4f,e=v=r; u=diag(N,A,A,A)r,N,A>0, with B and common endpoint arbitrary where such actual solutions exist",
            "exact_speed": "c_T=A/N=1+A n_h/(4q)>1 when n_h>0",
            "physical_TT_coefficients": "G_T,h=2GN/A³>0,F_T,h=2G/(NA)>0; F_T,h/G_T,h=A²/N²",
            "source_and_projection": "gamma=(gamma_e+gamma_v)/2,delta=(gamma_e-gamma_v)/2; gamma_u=gamma+A j/q; physical source j gamma and contact A j²/(2q)",
            "time_dependence": "Exact exchange symmetry decouples common and odd tensor channels despite time-dependent A,N; unequal Einstein coefficients do not inherit this",
            "normalization": "f_c=sqrt(a_h³G_T,h)/2,Y=f_c gamma,theta=D_T log f_c; canonical boundary -(theta Y²)'/2 and frequency c_T²kphys²-f_c''/f_c",
            "relative_mass_caution": "Relative algebraic mass 2qA²/(GN)>0 is not an adiabatic gap or canonical-frequency theorem; common cosmological channel is not a stationary mass pole",
            "matching_scope": "The common physical tensor already exceeds the actual matter cone; integrating out only the relative channel cannot fix this exact symmetric quadratic action",
        },
        "nonempty_actual_local_solution": {
            "same_specified_extension": "B=-6q,beta4g=beta4f=-q,G=F,q,G,epsilon>0; one free canonical scalar V=0",
            "domain_and_clock": "A=1+x>1,a_r>0,t common-r proper time; dT=N dt physical h time",
            "constraints": "rho=p=12q(A-1)/(epsilon A),N=A/(6A-5),H_r²=2q(A³-1)/(3G)",
            "analytic_ODE": "A'=-6H_r A(A-1)/(6A-5),a_r'=H_r a_r,psi'=N sqrt(2rho),H_r=positive square root",
            "full_equations": "H_r'=-6qA³(A-1)/[G(6A-5)]; both Einstein lapses/accelerations, all u components and actual free-scalar current equation vanish",
            "existence": "Analytic local ODE and open positive domain give a local actual rolling family for every A0>1,a_r0>0; no numeric uniform time radius computed",
            "initial_A2": {key: str(value) for key, value in rolling.initial_data().items()},
            "scope": "H_h=H_r/A>0, so this non-vacuous example is expanding, not CD, a bounce, global solution, scalar/vector-health or cutoff witness",
        },
        "quantitative_conditional_error_budget": {
            "exact_gap": "delta=A n_h/(4q); c_T²-1=delta(2+delta),F_T,h-G_T,h=G_T,h delta(2+delta)",
            "coefficient_threshold": "If |deltaG|<=E_G,|deltaF|<=E_F,G_T,h-E_G>0, then Fnew-Gnew>=G_T,h delta(2+delta)-E_G-E_F",
            "meaning": "A positive margin preserves the faster cone; an error sum at least the uncorrected gap is necessary, not sufficient, to repair it in that physical two-derivative representative",
            "numeric_interface": "Exact finite real inputs only; reject bool, floats, nonfinite values and unstated/wrong sign domains",
            "conditional_M1_comparison": "If the frozen canonical free-chi velocity M_target/(10tau) is preserved, other matter is positive canonical, and 2G=M_target² with m0²=2q/G, then c_T-1>=A/[100(m0 tau)²]",
            "not_an_error_estimate": "No actual omitted-operator, loop, background-displacement, nonlocal/dispersive or cutoff bound is supplied; one-scalar example is not identified with old clock-plus-chi matter",
        },
        "countercontrols_scope": {
            "mixed_sign_actual_flat": "p_g=-2,p_f=1,B=3,beta4g=2,beta4f=-1,e=v=u=I,zero matter solves all Euler equations; u=2e-v and literal flat relative spring=4>0",
            "mixed_sign_flat_limit": "Positive mass 4(1/G+1/F) and positive EH kinetics are flat quadratic controls only, not full rolling/quantum health",
            "mixed_sign_auxiliary_only": "n_e=1/2,n_v=n_u=a_e=a_v=a_u=1,p_g=-2,p_f=1,B=5/2,epsilon=1,rho=p=1 gives c_e=1/2,c_f=1; only u equations, e/v fail for these static frames",
            "nonrolling": "n_h=0 has no strict positive gap; the positive-link flat vacuum remains a valid control",
        },
        "primary_source": {"url": "https://arxiv.org/pdf/1804.04671", "used_for": "Action/source equations 2.1-2.17 and actual-vierbein warning section 5.1, not a rolling TT formula or imported cutoff"},
        "verification_boundary": "Exact symbolic and independent Fraction/polynomial replay plus separately authored covariant/action/source audit and written local-existence/positive-weight proof; not proof-assistant formalized",
        "not_established": [
            "A universal asymmetric light-only matching exclusion from a full-parent high-frequency cone",
            "The common-channel theorem on an exchange-breaking branch, mixed/negative links, singular lapses, or noncanonical/quantum anisotropic sources",
            "A stationary massive pole or positive heavy Green function inferred from a cosmological algebraic mass",
            "An original CD/M1 background/operator match, or scalar/vector stability of the auxiliary extension",
            "A bound on finite-band time advance, higher-operator or loop correction, cutoff or UV completion",
            "A closed causal curve or universal UV inconsistency solely from a cone comparison",
            "An exclusion across changed actions, additional u kinetic terms, or all trimetric/composite possibilities",
            "Closure of S6, either original P8 track, or P8",
        ],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Actual auxiliary-metric tensor-cone certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.14.AUXILIARY: actual physical common tensor-cone screen replay passed; general matching and P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
