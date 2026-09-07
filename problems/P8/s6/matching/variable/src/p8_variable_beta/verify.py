"""Read-only sourced variable-beta local solution and exact inner-TT replay."""

import argparse
import hashlib
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

import sympy as sp
from p8_hr_tree import verify as prior

from . import background, canonical, independent, inner, model, tensors

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"variable-beta-local.json"
P8 = next(path for path in ROOT.parents if path.name == "P8")
TREE_SHA = "d3d4c612fd0942bc66f0005ed21b76d23ca86285a461d87b32e1a2d79f242bb2"
CONTRACT_SHA = "d00df35d5821b05baa61e897725a22ed1ac21d0536b74bdaece3ab0093215901"
CONVENTION_SHA = "73dd18f7a4b56a0205f9fbc6a4b09b213a3c280b3f0274414a255635b04e2b09"
CD_SHA = "caf8c8e688a7565b9d00f921c099a28da00f97522ed26ad182a8227eb80cd4dd"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serialize(value):
    """Exact JSON data only; never turn a rounded number into a bound."""
    if isinstance(value, (Fraction, sp.Basic)):
        if isinstance(value, sp.Basic) and value.has(sp.Float, sp.oo, sp.zoo, sp.nan):
            raise TypeError("Only exact finite report expressions are allowed")
        return str(value)
    if isinstance(value, dict):
        return {str(key): serialize(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serialize(item) for item in value]
    if value is None or isinstance(value, (str, int, bool)):
        return value
    raise TypeError("Only exact report values are allowed")


@cache
def prior_checks():
    pins = {"S6_17_constant_beta_tree_scope_and_replay": (prior.REPORT, TREE_SHA),
            "adopted_S6_matching_contract": (P8/"s6"/"FORMULATION.md", CONTRACT_SHA),
            "P8_B_curvature_dictionary": (P8/"FORMULATION.md", CONVENTION_SHA),
            "original_CD_matter_trajectory_context": (P8/"certificates"/"witness-CD_matter.json", CD_SHA)}
    if any(sha(path) != expected for path, expected in pins.values()):
        raise ValueError("A frozen scope, physical convention, trajectory or contract changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {key: expected for key, (_, expected) in pins.items()}


def independent_bridges():
    """Compare independent arithmetic profiles and differentiated jets, not verdicts."""
    solutions = []
    for fixture in independent.background_fixtures():
        actual = background.evaluate(fixture["u"], fixture["c"])
        if any(actual[key] != value for key, value in fixture["values"].items()):
            raise ValueError("Primary background versus independent Fraction/Taylor construction failed")
        solutions.append({"u": str(fixture["u"]), "c": str(fixture["c"]),
                          "primary_minus_independent": dict.fromkeys(fixture["values"], "0")})
    jets = []
    d = canonical.derive()
    u, c = d["u"], d["c"]
    center = {key: sp.factor(d[key].subs(u, 0)) for key in
              ("K1", "K2", "mass_squared", "N_sum", "N_relative", "c_light_squared", "c_heavy_squared")}
    center["omega_prime_squared"] = sp.factor(sp.diff(d["omega"], u).subs(u, 0)**2)
    center["D_squared"] = sp.factor(d["D"].subs(u, 0)**2)
    center["mass_variation_ratio"] = sp.factor(sp.diff(d["mass_squared"], u, 2).subs(u, 0)/center["mass_squared"]**2)
    for fixture in independent.canonical_fixtures():
        expected = {key: value for key, value in fixture.items() if key != "c"}
        actual = {key: sp.cancel(value.subs(c, fixture["c"])) for key, value in center.items()}
        if actual != expected:
            raise ValueError("Independent true Taylor jets versus primary canonical derivatives failed")
        jets.append({"c": str(fixture["c"]), "primary_minus_independent": dict.fromkeys(expected, "0")})
    return {"actual_background_profiles": solutions, "canonical_second_jets": jets}


def checked_controls():
    bg, tt, forcing = background.controls(), tensors.controls(), inner.controls()
    if bg["omitted_clock_exchange_at_c1_u1over10"] == 0:
        raise ValueError("The old constant-coefficient Bianchi omission did not fire")
    if (bg["c1_center"]["kbar"], bg["c4_center"]["kbar"]) != (sp.Rational(5599, 100), sp.Rational(799, 100)):
        raise ValueError("The actual two canonical-clock center controls changed")
    if (tt["c1_center_shift_coefficient"], tt["c4_center_shift_coefficient"]) != (-sp.Rational(128, 3), sp.Rational(8, 3)):
        raise ValueError("The literal opposite vector-kinetic signs changed")
    if (tt["c1_center_f_tensor_speed_squared"], tt["c4_center_f_tensor_speed_squared"]) != (sp.Rational(1, 4), 4):
        raise ValueError("The actual physical-clock full tensor cones changed")
    if forcing["omitted_derivative_residual"] != sp.Rational(1, 5) or forcing["particular_to_algebraic_ratio"] != sp.Rational(5, 6):
        raise ValueError("The limiting, nonphysical-source algebraic omission control changed")
    d = canonical.derive()
    correction = inner.canonical_adiabatic_checks()
    if sp.simplify(correction["correction_at_center_limit"]-d["kbar"]**2+sp.Rational(18, 5)) != 0:
        raise ValueError("The exact canonical finite correction changed")
    if sp.simplify(correction["correction_second_at_center_limit"]+sp.Rational(24, 5)*d["kbar"]**2+sp.Rational(252, 5)) != 0:
        raise ValueError("The exact canonical second-derivative correction changed")
    invalid = [lambda value=value: background.evaluate(0, value)
               for value in (0, 2, sp.Rational(3, 2), 5, True, 0.1, sp.Float("0.1"), sp.oo, sp.nan)]
    invalid += [lambda value=value: background.evaluate(value, 3)
                for value in (sp.Rational(-11, 100), sp.Rational(11, 100), True, 0.1, sp.oo)]
    invalid += [lambda: independent.Jet(True), lambda: independent.Jet(1.0),
                lambda: independent.Jet(1)**True, lambda: independent.Jet(1)**-1,
                lambda: serialize(sp.Float("0.1")), lambda: serialize(sp.oo)]
    rejected = 0
    for call in invalid:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(invalid):
        raise ValueError("An invalid or inexact local-family input was admitted")
    return serialize({"source_aware_background_controls": bg, "literal_TT_shift_controls": tt,
                      "finite_canonical_correction_limits": correction,
                      "limiting_operator_omission_control": forcing,
                      "invalid_exact_domain_calls_rejected": rejected})


@cache
def build_report():
    previous = prior_checks()
    residuals = {"literal_source_aware_two_lapse_action": model.checks(),
                 "canonical_and_optional_conformal_source_bookkeeping": model.source_checks(),
                 "actual_metric_and_canonical_scalar_solution": background.checks(),
                 "literal_TT_elementary_potential_and_shift": tensors.checks(),
                 "literal_spatial_curvature_gradient": tensors.spatial_curvature_check(),
                 "exact_transverse_shift_elimination": tensors.vector_checks(),
                 "full_moving_weight_normalization": canonical.checks(),
                 "canonical_action_boundary_and_adjoint": canonical.action_checks(),
                 "actual_canonical_center_data": canonical.center_checks(),
                 "exact_inner_equation_and_nonadiabatic_limits": inner.checks()}
    if any(sp.simplify(value) != 0 for group in residuals.values() for value in group.values()):
        raise ValueError("A source, actual-equation, tensor or inner-operator identity failed")
    margins = background.domain_checks()
    reserved = ROOT/"tests"/"test_variable_independent_audit.py"
    if not reserved.is_file():
        raise ValueError("The separately authored covariant audit must exist before freezing")
    sources = sorted(ROOT.glob("src/p8_variable_beta/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1, "claim": "P8-S6.20.VARIABLE", "date": "2026-09-06",
        "status": "ACTUAL_LOCAL_VARIABLE_BETA_CD_TRAJECTORY_AND_NONADIABATIC_INNER_TENSOR; ORIGINAL_P8_OPEN",
        "prior_context_sha256": previous,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
        "formulation": "FORMULATION.md", "written_proof": "notes/proof.md", "source_audit": "notes/sources.md",
        "exact_residuals": {name: dict.fromkeys(group, "0") for name, group in residuals.items()},
        "strict_rational_clock_domain_margins": serialize(margins),
        "independent_Fraction_replay": serialize(independent.checks()),
        "primary_independent_bridges": independent_bridges(), "checked_controls": checked_controls(),
        "action_and_physical_frame": {
            "Einstein_terms": "-M²/2 integral(sqrt|g| R_B[g]+sqrt|f| R_B[f]); both constant positive Einstein coefficients M²",
            "signature_curvature": "+---; R_B=-6(DH+2H²); 3M² H²=rho; -2M² DH=rho+p",
            "actual_matter_metric": "g, with actual proper clock dT=N_g dt; no implicit conformal or composite frame",
            "two_literal_scalar_actions": "+1/2 integral sqrt|g|[(partial phi)²_g+(partial chi)²_g]; phi interaction-sourced, chi exactly free",
            "interaction": "-2 integral sqrt|g| sum beta_n(phi) e_n(sqrt(g^-1 f)); positive diagonal FLRW root",
            "beta_units": "beta_n(phi)=M²/tau²*b_n(u(phi/M),c); functions of the dynamical clock, never external beta(t)",
            "optional_f_scalar": "An absent or constant free f scalar contributes no source; beta is independent of it",
        },
        "undivided_source_identities": {
            "U_V_P": "U=beta0+3beta1*y+3beta2*y²+beta3*y³; V=beta1+3beta2*y+3beta3*y²+beta4*y³; P=2(beta1+2beta2*y+beta3*y²)",
            "literal_nulls": "nu_Ig=(y-c)P; nu_If=(c-y)P/(c*y³)",
            "clock_forces": "F_g=-2(U_phi_g+c V_phi_g); F_f=-2(U_phi_f+c V_phi_f)/(c*y³)",
            "combined_balance": "C_g=c B; C_f=-B/(c*y³); B=3P(yH_f-H_g)+2(U_phi_f*q_f-V_phi_g*q_g)",
            "consequence": "On-shell B=0, not the old velocity factor. No division by P,H,or clock speed in the generic identity",
            "optional_conformal_dictionary": "For h=A²g, matter balance is alpha*q*(rho_E-3p_E), clock force its negative; physical H_h=A^-1(H_g+alpha*q); the actual fixture uses A=1",
        },
        "actual_local_solution": {
            "domain": "M,tau>0; |u|<=1/10, u=T/tau; c=1 or any fixed 2<c<=4; no uniform bounded beta as c->2+",
            "geometry": "N_g=1,N_f=c,a=(1+u²)²,b=2/(1+u²)²,y=2/(1+u²)^4,H_g=h/tau,H_f=-h/(c*tau),h=4u/(1+u²)",
            "positive_canonical_clocks": "nbar=2(y³/c-1)h'; chi_T=M/(10tau*a³); kbar=nbar-1/(100a^6); phi=M integral_0^u sqrt(kbar)dv",
            "profiles": "b1=y³h'/[c(c-y)]; b2=b3=0; b4=3h²/(2c²)-b1/y³; b0=(3h²-nbar/2)/2-3b1*y",
            "equations": "Both lapse equations, both null/spatial equations, phi and chi equations, and sourced undivided Bianchi hold exactly; remaining FLRW components vanish by symmetry",
            "analytic_action": "Strict positive kbar plus denominator separation gives an analytic inverse clock and analytic beta_n(phi) on each fixed-c local range",
            "uniform_clock_bounds": "kbar>3599/100 for c=1; kbar>449/100 for 2<c<=4; bare clock positivity is not reduced coupled-scalar health",
            "original_target_boundary": "Locally the same CD scale factor and free-chi trajectory, not the original DHOST clock/action/physical operator matching",
        },
        "literal_tensor_and_vector_result": {
            "unit_TT_polarization": "e_ij e_ij=1; S_T=1/8 integral dT a³{M²[gamma_gdot²-q²gamma_g²]+M²y³/c[gamma_fdot²-(c/y)²q²gamma_f²]-mu(gamma_f-gamma_g)²}",
            "mu": "mu=2y[beta1+beta2(c+y)+beta3*c*y]; mu=yP for this beta1-only family",
            "scalar_TT_response": "Homogeneous canonical scalar action has unit determinant under traceless exponential spatial perturbations; no independent TT matter response assumed",
            "literal_relative_shift_coefficient": "a^5*y²*P/[2N_g(c+y)]",
            "vector_schur": "K_relative=AB*C_s*k²/[AB*k²+(A+B)C_s], A,B>0; negative at sufficiently high formal k for c=1; positive kinetic for c>2, not a full vector-health claim",
            "full_principal_cones": "Squared speeds 1 and c²/y² relative to actual g; second is wider on 2<c<=4",
            "locked_coefficient_only": "c_L²=(1+c*y)/(1+y³/c), not a controlled integrated-out physical light mode",
            "algebraic_mass": "tau² m_alg²=2y*h'*(y³+c)/[c(c-y)]; neither its sign nor a center canonical diagonal is a rolling gap theorem",
        },
        "complete_canonical_map": {
            "weights_M_equals_one": "K1=a³/8,K2=b³/(8c),Ksum=K1+K2,KR=K1K2/Ksum,w_i=K_i/Ksum; fsum=sqrt(2Ksum),fR=sqrt(2KR)",
            "coordinates": "l=fsum*(w1 gamma_g+w2 gamma_f),Q=fR*(gamma_f-gamma_g); restore both f by M; action in u has overall 1/tau",
            "actual_inverse": "gamma_g=l/fsum-w2*Q/fR; gamma_f=l/fsum+w1*Q/fR; retain this map for physical data/source comparisons",
            "connections": "theta_i=f_i'/f_i,N_i=theta_i'+theta_i²,omega=w2'/(2sqrt(w1w2)); all primes d/du",
            "boundary": "-(theta_sum*l²+theta_relative*Q²)'/2",
            "diagonals": "V_LL=qbar²*c_L²-N_sum; V_HH=mbar²+qbar²*c_H²-N_relative-4omega²",
            "mixed_operators": "B=2omega partial_u-2omega theta_sum+qbar²D; B*=-2omega partial_u-2omega'-2omega theta_sum+qbar²D",
            "future_retarded_data": "A retarded inverse would give Q=Q_hom-G_R B l; no homogeneous-data or physical-source suppression supplied here",
        },
        "shrinking_inner_window_theorem": {
            "scaling": "c=2+epsilon²,u=epsilon*x,epsilon>0; fixed finite x interval and bounded fixed dimensionless momenta kappa=tau*kcom",
            "exact_mass_profile": "mbar²=16(1-u²)(8+c*d^12)/[c*d^14*(c*d^4-2)]",
            "center_mass": "mbar0²=16(c+8)/[c(c-2)] ~80/(c-2)",
            "exact_second_ratio": "(mbar²)''0/(mbar0²)²=-c(7c²+146c-240)/[8(c+8)²] ->-1/5",
            "full_canonical_ratios": "V_HH-mbar² and its second derivative are analytic finite at (u,c)=(0,2), so full V_HH has the same -1/5 ratio and positive-square-root -1/10 ratio for small positive epsilon",
            "limiting_equations": "l_xx=0; Q_xx+80/(1+8x²)Q=0; all scaled off-diagonal potential and derivative coefficients tend to zero",
            "compact_continuation": "After exact epsilon² denominator removal all scaled coefficients extend smoothly to epsilon0 on fixed compact x sets; uniform first-order coefficient convergence gives Cauchy-propagator convergence for convergent data",
            "explicit_local_domain_choice": "0<epsilon<=min(1,1/[10 max(1,X)]) keeps |x|<=X inside the certified actual-clock interval",
            "physical_projection_limit": "gamma_g ->2(l-2Q)/(M sqrt5); nonzero relative projection with finite canonical normalization, not a source-response assertion",
            "limiting_action_guard": "epsilon0 is not a regular finite-coefficient parent; beta profiles diverge",
            "window_guard": "Only |T|<=epsilon*tau*X. No computed epsilon threshold, fixed physical/CD window, growing momentum band or retarded physical matching error",
            "forcing_control": "For limiting constant normalized forcing, R=(1+8x²)/96 is exact; algebraic R0=(1+8x²)/80 has residual1/5 and particular ratio5/6. No equal retarded data or physical TT source was imposed",
        },
        "primary_source": "https://arxiv.org/pdf/2501.16442; section4.1 and AppendixB, with source/frame/interaction dictionary in notes/sources.md; no novelty, cutoff or model-realization claim",
        "verification_boundary": "Exact identities, independent coefficient/Fraction/Taylor outputs and separately authored audits support the written inverse-function and compact-ODE proof; no proof-assistant formalization",
        "not_established": [
            "No full reduced scalar constraints, scalar/vector health, or physical subcutoff growth-rate conclusion",
            "No rolling spectral gap or adiabatic elimination inferred from algebraic mass height or a frozen canonical diagonal",
            "No light-only matching exclusion from the wider full tensor cone, locked coefficient, or limiting normalized forcing",
            "No complete-CD solution or uniform coefficient/physical-window limit as c->2+",
            "No original DHOST C/D operator and physical-source/data map; no loop, cutoff, positivity or UV-completion result",
            "No change to completed scoped photon objective or frozen linear classification; original S6, P8(b) and P8 remain open",
        ],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Variable-beta local/tensor certificate differs from exact read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.20.VARIABLE: actual local background and nonadiabatic inner tensor replay passed; original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
