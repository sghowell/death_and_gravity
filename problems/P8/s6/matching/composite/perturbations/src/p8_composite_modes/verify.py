"""Read-only source-hashed tensor/vector screen and immutable lineage replay."""

import argparse
import hashlib
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

import sympy as sp
from p8_composite import verify as prior

from . import independent, jets, model, shift, tensor, vector

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"composite-modes.json"
PRIOR_SHA = "7708eebb7f10be8bad523029aa7a53921c943fd6f7df66691b71651d090ec06f"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@cache
def prior_checks():
    if sha(prior.REPORT) != PRIOR_SHA:
        raise ValueError("Pinned S6.6 composite-background certificate changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return PRIOR_SHA


def exact_json(value):
    if isinstance(value, dict):
        return {key: exact_json(item) for key, item in value.items()}
    if isinstance(value, (Fraction, sp.Basic)):
        return str(value)
    return value


def residuals():
    groups = {"literal_source_aware_TT": tensor.checks(),
              "literal_two_by_two_shift_and_effective_metric": shift.checks(),
              "full_vector_constraints_and_canonical_time_boundary": vector.checks(),
              "actual_clock_generic_acceleration_jets": jets.checks(),
              "independent_physical_scale_restoration": jets.scale_checks()}
    for name, values in groups.items():
        if any(sp.simplify(value) != 0 for value in values.values()):
            raise ValueError("A composite-mode identity failed: "+name)
    return groups


@cache
def build_report():
    previous = prior_checks()
    groups = residuals()
    rational = independent.checks()
    positives = shift.positive_factors()
    vd = vector.derive()
    positives.update({"positive_K_at_mu_zero": vd["K"].subs(vd["mu"], 0),
                      "positive_physical_K_at_mu_zero": vd["K_T"].subs(vd["mu"], 0)})
    if any(value.is_positive is not True for value in positives.values()):
        raise ValueError("A stipulated-domain positive coefficient failed")
    controls = {module.__name__.rsplit(".", 1)[-1]: module.negative_controls()
                for module in (tensor, shift, vector)}
    if any(sp.simplify(value) == 0 for values in controls.values() for value in values.values()):
        raise ValueError("A genuine omission control failed to fire")
    signs = jets.negative_controls()
    if not (signs["free_frozen_negative_speed_coefficient"] < 0
            and signs["CD_frozen_negative_speed_coefficient"] < 0
            and signs["independent_scale_reverses_leading_sign"] > 0
            and signs["threshold_does_not_determine_higher_jets"] == 0):
        raise ValueError("A frozen-candidate exclusion or independent-scale guard failed")
    source = vector.source_normalization_checks()
    if not (source["4_11_C_matches_literal"] == 0
            and source["printed_4_12_K_over_literal"] == 2
            and sp.cancel(source["printed_4_12_U_over_literal"]
                          -2*model.F*model.y**2/(model.G+model.F*model.y**2)) == 0
            and sp.cancel(source["printed_frequency_over_literal"]
                          -model.F*model.y**2/(model.G+model.F*model.y**2)) == 0):
        raise ValueError("Literal 4.8 source-normalization audit failed")
    jd = jets.derive()
    sources = sorted(ROOT.glob("src/p8_composite_modes/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {"schema": 1, "claim": "P8-S6.7.COMPOSITE", "date": "2026-09-06",
            "status": "TWO_FROZEN_COMPOSITE_BOUNCES_FAIL_LOCAL_FORMAL_VECTOR_PRINCIPAL_HEALTH; INDEPENDENT_SCALE_CAN_REVERSE_LEADING_SIGN; GLOBAL_SCALAR_CUTOFF_MATCHING_OPEN",
            "prior_S6_6_composite_background_sha256": previous,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
            "formulation": "FORMULATION.md", "written_proof": "notes/modes.md",
            "primary_source_audit": "notes/sources.md",
            "primary_sources": {
                "unreduced_TT_vector_action_and_dictionary": "https://arxiv.org/pdf/1501.02790",
                "composite_matter_coupling": "https://arxiv.org/pdf/1408.1678",
                "generic_constraint_and_EFT_caution": "https://arxiv.org/pdf/1409.3834",
                "specified_flat_normalization_cutoff_upper_bound": "https://arxiv.org/pdf/1506.00666",
                "higher_source_order_completion": "https://arxiv.org/pdf/1804.04671"},
            "exact_residuals": {name: dict.fromkeys(values, "0") for name, values in groups.items()},
            "independent_Fraction_coefficients_and_controls": rational,
            "positive_domain_coefficients": exact_json(positives),
            "nonzero_omission_controls": exact_json(controls),
            "source_normalization_audit": exact_json(source),
            "exact_jet_and_threshold_controls": exact_json(signs),
            "domain": {
                "action": "Pinned HR plus one canonical scalar on g_eff=g(alpha I+beta sqrt(g^-1 f))^2; +---; G,F,m4,alpha,beta>0 and constant real beta_n",
                "background": "Full source-aware background equations; common flat positive-root FLRW, smooth positive finite a,b,Ng,Nf",
                "physical_frame": "Ae=a*(alpha+beta*y),Ne=Ng*(alpha+beta*c),dT=Ne*dt; physical metric is g_eff, not old CD/M1 g",
                "source": "rho+p=(D_eff chi)^2; induced g/f sources are not separately conserved",
                "vector_chart": "k>0 fixed nonzero comoving wavenumber, Xi>0; sufficiently small transverse shifts on the positive analytic square-root branch",
                "sufficient_Xi_domain": "P>0,rho>0,rho+p>0; all hold near the two frozen bounce points; no arbitrary-beta positivity theorem"},
            "TT": {"mu": str(tensor.derive()["mu"]),
                   "normalization": "TT contraction prefactor1/8; literal diagonal exp D=diag(1,-1,0) has norm2 and relative L2=-Ng*a^3*mu*(gamma-h)^2/4",
                   "kinetic": "G*a^3/(8Ng),F*b^3/(8Nf)>0 per TT contraction",
                   "mass_boundary": "mu is algebraic relative stiffness, not a computed rolling canonical gap or growth rate"},
            "finite_shift": {"Z_v": str(shift.derive()["Z_v"]), "Xi": str(shift.derive()["Xi"]),
                             "root": "S2=(M+c*y I)/sqrt[(c+y)^2-w],w=b^2(sigma_f-sigma_g)^2/Ng^2",
                             "matter_rule": "Hold chi_dot fixed; vary both sqrt|g_eff| and g_eff^00",
                             "zero_mu": "Xi is evaluated by its algebraic expression, never by dividing mu/cV^2 at mu=0"},
            "vector": {"field": "q=k(E-S),spatial vector convention h_ij=i(k_i E_j+k_j E_i)/2",
                       "completed_square": "A0=k^2/2,B0=F*y^3*k^2/(2G*c),C0=a^2*Xi/(2G); harmonic inertia=(A0^-1+B0^-1+C0^-1)^-1",
                       "C": str(vd["C"]), "K": str(vd["K"]), "U": str(vd["U"]),
                       "K_T": str(vd["K_T"]), "U_T": str(vd["U_T"]),
                       "physical_principal_speed_squared": str(vd["physical_principal_speed_squared"]),
                       "canonical_normalization": "Q=sqrt(2K_T)q,Omega_T^2=U_T/K_T-(sqrt(K_T))_TT/sqrt(K_T); full time boundary retained",
                       "mu_zero_boundary": "K_T remains positive; mu=0 alone establishes neither vanishing physical gap nor strong coupling"},
            "analytic_bounce_screen": {
                "parameters": "G=F=M^2,m4=M^2*m^2,alpha=beta=1,beta_n=(0,0,1,0,0); u=m*T; rho,p divided by M^2*m^2",
                "initial": "y=c=1,rho_bar=p_bar=1/2,Ng=Nf=1/2; X=-Y=sqrt(7/3); A=hprime(0) with h=H_eff/m",
                "Xprime": str(jd["X_prime"]), "Yprime": str(jd["Y_prime"]),
                "yprime": str(jd["y_prime"]), "cprime": str(jd["c_prime"]),
                "mu_over_m4_u2": str(jd["mu_quadratic"]), "Xi_over_m4_at_zero": str(jd["Xi_initial"]),
                "physical_speed_u2": str(jd["physical_speed_quadratic"]),
                "free": "A=7/36,mu/m4=-u^2/18+O(u^3),cVeff^2=-u^2/108+O(u^3)",
                "frozen_CD": "A=4,mu/m4=-23u^2/3+O(u^3),cVeff^2=-23u^2/18+O(u^3)",
                "proof": "Analyticity and the exact nonzero negative quadratic coefficient imply a two-sided punctured neighborhood with mu<0,cVeff^2<0 while Xi,K_T remain positive",
                "independent_scale_guard": "For CD h=4u/[(m*tau)^2+u^2],A=4/(m*tau)^2 and mu/m4=[(m*tau)^2-24]u^2/[3(m*tau)^2]+O(u^3). m*tau>sqrt24 reverses only the leading sign; equality requires higher jets"},
            "verification_boundary": [
                "Literal TT source variation and 3D spatial curvature; full 2x2 root/effective-metric source variation",
                "Exact finite-k two-shift constraint elimination, not a singular mu/cV chart at the bounce",
                "Physical clock and ruler, k-dependent inertia and normalization time boundary retained",
                "Independent Fraction coefficient normal forms, exact sign controls, root-owned independent covariant tests",
                "Pinned background action and full prior lineage replayed without modifications"],
            "not_established": [
                "A numerical width of the negative-principal interval, a quantitative accumulated growth bound or large instability below a known cutoff",
                "Full coupled scalar constraints, kinetic/gradient matrices, finite-band evolution, nonlinear interactions or a background-dependent BD-mode cutoff",
                "A physical rolling tensor/vector heavy-mode gap, or vanishing such a gap merely from mu=0",
                "Health for independently rescaled CD beyond the sign of the leading vector coefficient; the threshold equality is undecided",
                "A globally controlled CD solution, old free-M1 matter-frame matching, positivity closure or S6/P8 completion",
                "A quantum front-velocity verdict or a general exclusion of composite, bimetric or ultraviolet-completed theories"]}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Composite-mode certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.7.COMPOSITE: exact tensor/vector constraints and frozen local vector-principal obstruction replay passed; general scalar/cutoff matching OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
