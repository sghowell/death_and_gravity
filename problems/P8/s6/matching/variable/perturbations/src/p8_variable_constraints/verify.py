"""Read-only regular scalar constraints/vector action replay; scalar health OPEN."""

import argparse
import hashlib
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

import sympy as sp
from p8_variable_beta import verify as prior

from . import action, auxiliary, domain, independent, observables, reduction, vector

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"variable-constraints.json"
PRIOR_SHA = "335cd52028baf56b30c75377aa7e6edd7ec86db2799f49b13f467242674fc4d0"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serialize(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, sp.Basic):
        if value.has(sp.Float, sp.oo, -sp.oo, sp.zoo, sp.nan):
            raise TypeError("Only exact finite report expressions are permitted")
        return str(value)
    if isinstance(value, dict):
        return {str(key): serialize(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serialize(item) for item in value]
    if value is None or isinstance(value, (str, int, bool)):
        return value
    raise TypeError("Only exact JSON/report values are permitted")


@cache
def prior_checks():
    if sha(prior.REPORT) != PRIOR_SHA:
        raise ValueError("The frozen S6.20 actual action/background certificate changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {"S6_20_actual_variable_beta_action_background_and_adopted_contract": PRIOR_SHA}


@cache
def independent_bridges():
    phase = []
    for fixture in independent.scalar_fixtures():
        d = reduction.jets(fixture["u"], fixture["c"], 1)
        at = {action.K: fixture["K"]}
        count = 0
        for name in ("M", "H"):
            for order in range(2):
                actual = d[name][order].subs(at).applyfunc(sp.cancel)
                expected = sp.Matrix(fixture[f"{name}{order}"])
                if actual != expected:
                    raise ValueError("Independent Fraction polarization disagreed with the full phase Hessian jet")
                count += actual.rows*actual.cols
        for order in range(2):
            if d["p0"][order].subs(at).applyfunc(sp.cancel) != sp.Matrix([fixture["p0"][order]]):
                raise ValueError("Independent secondary momentum jet failed")
            if sp.cancel(d["D"][order].subs(at)) != fixture["D"][order]:
                raise ValueError("Independent secondary denominator jet failed")
            count += 7
        phase.append({**{key: str(fixture[key]) for key in ("u", "c", "K")},
                      "exact_scalar_coefficient_comparisons": count, "all_residuals_zero": True})
    vd = vector.derive()
    vectors = []
    for fixture in independent.vector_fixtures():
        at = {vd["u"]: fixture["u"], vd["c"]: fixture["c"], action.K: fixture["K"]}
        for independent_key, primary_key in (("kinetic", "kinetic"), ("bare", "omega_bare_squared"),
                                             ("pump", "normalization_pump"), ("speed_squared", "speed_squared"),
                                             ("mass_squared", "mass_algebraic_squared")):
            if sp.cancel(vd[primary_key].subs(at)) != fixture[independent_key]:
                raise ValueError("Independent literal vector-shift/Taylor output failed")
        vectors.append({**{key: str(fixture[key]) for key in ("u", "c", "K")},
                        "exact_scalar_coefficient_comparisons": 5, "all_residuals_zero": True})
    own, other = domain.derive(), independent.margin_coefficients()
    if own["c1_coefficients"] != other["c1"] or own["positive_branch_coefficients"] != other["positive"]:
        raise ValueError("Independent cleared-polynomial and Bernstein construction disagreed")
    coefficient_count = sum(len(values) for values in other.values())
    return {"full_regular_phase_jets": phase, "literal_vector_shift_and_pump": vectors,
            "independent_continuous_domain_coefficients": coefficient_count,
            "total_exact_output_comparisons": sum(item["exact_scalar_coefficient_comparisons"] for item in phase+vectors)
                                                +coefficient_count}


def checked_controls():
    toy = observables.nonuniform_controls()
    if toy["exact_transformed_free_oscillator"] != 0:
        raise ValueError("The nonuniform-map countercontrol failed")
    if toy["original_frozen_frequency_squared"] != -toy["transformed_center_frozen_frequency_squared"]:
        raise ValueError("The false frozen-frequency inference was not rejected")
    r, d = reduction.derive(), action.derive()
    at = {d["u"]: 0, d["c"]: 4, action.K: 1}
    correction = sp.diff(r["H_time_correction"], reduction.P0, 2).subs(at)
    actual = reduction.jets(0, 4, 0)["D"][0]
    if correction != sp.Rational(7, 144) or actual != -sp.Rational(5, 24):
        raise ValueError("The canonical time-boundary omission control changed")
    # This negative block is an explicit guard against calling a canonical
    # instantaneous energy block the physical kinetic/signature theorem.
    old_minor = sp.factor(reduction.jets(0, 4, 0)["H"][0][3:5, 3:5].det()).subs(action.K, 1)
    if old_minor != -sp.Rational(17, 384):
        raise ValueError("The chart-dependent momentum-block control changed")
    oj = observables.observable_jets(sp.Rational(1, 100), 4)["O"][0]
    bracket = (oj*reduction.J*oj.T).applyfunc(sp.cancel)
    if bracket[0, 2] == 0:
        raise ValueError("The off-center physical-observable noncommutativity control vanished")
    # Warm exact caches before invalid inputs: float/bool aliasing must not
    # bypass an input boundary after a valid call.
    observables.equation(0, 4)
    auxiliary.derive(0, 4)
    invalid = [lambda args=args: reduction.jets(*args) for args in
               ((True, 4, 0), (0.0, 4, 0), (0, sp.Float(4), 0), (0, sp.oo, 0),
                (0, 2, 0), (sp.Rational(11, 100), 4, 0), (0, 4, True),
                (0, 4, 1.0), (0, 4, 3), (None, 4, 0), (0, sp.I, 0))]
    invalid += [lambda: observables.equation(0.0, 4), lambda: auxiliary.derive(0.0, 4),
                lambda: independent.scalar_jet_fixture(True, 4, 1),
                lambda: independent.scalar_jet_fixture(0, 4, 0),
                lambda: independent.inverse([[1, 1], [1, 1]]),
                lambda: serialize(sp.Float("0.1")), lambda: serialize(sp.oo)]
    rejected = 0
    for call in invalid:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(invalid):
        raise ValueError("An inexact, singular or out-of-domain input bypassed validation")
    return {"exact_nonuniform_free_oscillator": serialize(toy),
            "omitted_time_boundary_D_minus_actual_at_c4_K1": str(-correction),
            "instantaneous_canonical_two_momentum_minor_NOT_ghost_verdict": str(old_minor),
            "off_center_physical_q_Dirac_bracket_xi_Psi_f": str(bracket[0, 2]),
            "invalid_or_inexact_calls_rejected_after_cache_warmup": rejected}


@cache
def build_report():
    previous = prior_checks()
    residuals = {"literal_regular_canonical_constraint_and_boundary": reduction.checks(),
                 "literal_transverse_geometry_potential_and_Shift_Schur": vector.checks(),
                 "actual_center_vector_clock_and_sign_controls": vector.center_checks(),
                 "fixed_K_compact_inner_vector_limit": vector.inner_checks(),
                 "finite_K_physical_center_and_nonuniform_map_controls": observables.center_checks(),
                 "continuous_domain_clock_and_margin_bridge": {
                     key: domain.derive()[key] for key in ("clock_F_identity", "margin_identity")}}
    for name, time, lapse in (("center_c4", 0, 4), ("future_c4", sp.Rational(1, 100), 4),
                              ("past_c3", sp.Rational(-1, 100), 3), ("center_c1_control", 0, 1)):
        residuals[f"all_unreduced_equations_{name}"] = auxiliary.checks(time, lapse)
    if any(sp.cancel(value) != 0 for group in residuals.values() for value in group.values()):
        raise ValueError("A literal action/constraint, physical-map, vector or equation identity failed")
    if not (ROOT/"tests"/"test_constraints_independent_audit.py").is_file():
        raise ValueError("The separately authored independent audit is required before freezing")
    sources = sorted(ROOT.glob("src/p8_variable_constraints/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    jets = observables.center_time_map()
    return {
        "schema": 1, "claim": "P8-S6.22.CONSTRAINTS", "date": "2026-09-06",
        "status": "EXACT_REGULAR_SCALAR_CONSTRAINTS_AND_VECTOR_DYNAMICS; SCALAR_HEALTH_AND_ORIGINAL_P8_OPEN",
        "prior_context_sha256": previous,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
        "formulation": "FORMULATION.md", "written_proof": "notes/proof.md", "source_audit": "notes/sources.md",
        "exact_residuals": {name: dict.fromkeys(group, "0") for name, group in residuals.items()},
        "exact_residual_count": sum(map(len, residuals.values())),
        "strict_continuous_secondary_constraint_margin": serialize(domain.checks()),
        "independent_Fraction_Taylor_polynomial_replay": serialize(independent.checks()),
        "primary_independent_output_bridges": independent_bridges(),
        "checked_omissions_and_invalid_inputs": checked_controls(),
        "literal_action_and_domain": {
            "physical_metric": "The same g and canonical phi/free chi of frozen S6.20; no external beta(t), conformal source, or G1 transfer",
            "curvature_EH": "+---, R_B=-6(DH+2H²), G_B00=3H²; Einstein action -M² R_B/2 on both metrics",
            "domain": "|u|<=1/10, c=1 or 2<c<=4 constant, M,tau>0, K=(tau*kcom)²>0; beta2=beta3=0",
            "unit_restoration": "u=T/tau,xbar=x/tau,phi_bar=phi/M,chi_bar=chi/M; overall action(Mtau)², physical frequency=omega_u/tau, spatial momentum=sqrtK/(tau*a)",
            "regular_gauge": "delta_phi=0 and E_g=0; requires positive actual clock and nonzero mode, never divides H",
            "EH_boundary": "-[a³h Sg²-b³h Sf²/c]'; all background derivatives retained before evaluating H=0",
            "full_action": "action.lagrangian retains both lapse and both shift equations; root literal ADM/square-root/Legendre audit and full unreduced-equation replay",
        },
        "regular_scalar_reduction": {
            "primary": "Cf=-h pf-2cbK qf+ca³P(Sg-Sf)=3ca³P Q",
            "canonical_map": "r=1+2bK/(3a³P),alpha=h/(3ca³P); qg=Q+rR+e/3+alpha PR,qf=R+alpha p0,pg=p0,pf=PR-rp0,pe=PE-p0/3",
            "boundary": "F=alpha p0 PR-r alpha p0²/2; add G=-r' p0R-alpha' p0PR+(r alpha'-alpha r')p0²/2 to H",
            "secondary_D": "-(1+c/y³)/(6a³)+h²(1+y/c)²/(a³k)-2h²y/(3c²a³P)+alpha'; exactly K independent",
            "strict_rank": "D<-1/8 throughout the stated box; D0=-5/24; one second-class pair leaves three canonical scalar pairs",
            "continuity": "For each fixed c and finite compact K band, regular canonical coefficients are continuous across the full local interval; no uniform scalar microlocal energy estimate inferred",
            "auxiliary_recovery": "ng=Sg-(A0-wpx)/(a³k),Bg=(3pEg-pg)/(2a³K),Bf=Bg-pe/(2Cs); nf from p0 preservation divided only by3ca³P",
        },
        "actual_physical_scalar_dictionary": {
            "coordinates": "xi=3(a³wz-PE)/(2aK),delta_phi_B=sqrt(k)*xi,chi_B=z+wxi,Psi_f_B=R+alpha p0-hxi",
            "Cauchy_map": "q=OZ; C=stack(O,O'+OA); q''=[O''+2O'A+O(A'+A²)] C^-1(q,q')",
            "conserved_form": "Z1^T Jcan Z2; pullback C^-T Jcan C^-1. This is not the instantaneous Hamiltonian momentum block",
            "center_invertibility": "F_K=c²(c-2)K²+4c(c-2)(c-8)K+1536(8-c)>0; detC=5F_K/[12K²(c-2)(6400-801c)] on2<c<=4,K>0",
            "finite_band_neighborhood": "Existence follows by compact continuity for each fixed finite positive K band and c; no momentum-uniform neighborhood or quantitative bound supplied",
        },
        "nonuniform_scalar_guard": {
            "c4_determinant_jets": serialize({key: jets[key] for key in ("det0", "det1", "det2", "relative_u2_K_coefficient")}),
            "inference_rejected": "Apparent center frozen polynomial and positive q-velocity bracket are NOT certified physical scalar cones/no-ghost results",
            "map_order_issue": "relative u²K coefficient5/6, off-center physical q Dirac noncommutativity, and nonuniform oscillator control require a uniform weighted scalar reduction",
            "next_mathematical_obligation": "Uniform scalar high-frequency kinetic/gradient/characteristic analysis with physical source/data maps; research work, not a user-intervention blocker",
        },
        "complete_vector_result": {
            "literal_transverse_action": "A K(Eg'-sg)²+B K(Ef'-sf)²+Cs(sf-sg)²-Dv K(Ef-Eg)²; A=a³/4,B=b³/(4c),Dv=a³mu/4,mu=yP",
            "exact_shift_elimination": "Kv=AB Cs K/[ABK+(A+B)Cs]",
            "physical_principal_speed": "cV²=(c+y)/(2y)>1 on2<c<=4; Kv>0 there for allK>0. Formal full-parent statement only",
            "negative_control": "c1 center Kv=-128K/[3(K-192)]; negative forK>192, singular atK192; not a below-cutoff growth or light-only exclusion",
            "normalization": "Positive branch V=sqrt(2Kv)(Ef-Eg),theta=Kv'/(2Kv); boundary -(theta V²)'/2; omega_can²=m_alg²+cV²K/a²-theta'-theta²",
            "center_c4": "Kv=8K/[3(K+16)],omega_can²=(9K²+382K+2112)/[6(K+16)]>0; not a rolling gap theorem",
            "inner": "For c=2+epsilon²,u=epsilon x,fixedK>0 and compactx: Kv->K/5,Vxx+80V/(1+8x²)=0; only shrinking physical window",
        },
        "verification_boundary": "Exact algebra, independently polarized Fraction/Taylor outputs, all140 independently derived Bernstein coefficients and separate source-aware ADM audits; continuity proof written, not proof-assistant formalized",
        "not_established": [
            "No scalar ghost/gradient/cone verdict from an instantaneous Hamiltonian or a nonuniform frozen observable polynomial",
            "No uniform scalar high-frequency energy estimate, finite-band scalar source matching or complete coupled health",
            "No rolling spectral gap, global CD parent solution, nonlinear stability, cutoff, loops, positivity or UV completion",
            "No original DHOST C/D operator match; no physical-frame or free-chi action change",
            "No dependence on separate S6.21 response or exploratory G1 action; no transfer between different lapse/clock jets",
            "Completed scoped photon objective and frozen linear classification unchanged; original P8 remains open",
        ],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Regular scalar/vector certificate differs from exact read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.22.CONSTRAINTS: regular scalar constraint and vector replay passed; scalar health and original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
