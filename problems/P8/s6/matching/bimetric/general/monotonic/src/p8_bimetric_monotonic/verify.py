"""Read-only normalized-Hubble monotonicity and sharp CD endpoint replay."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_bimetric_general import verify as prior

from . import controls, flow, fraction_checks

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"regular-hr-monotonic.json"
PRIOR_SHA = "e137b2392930eb9cec1abbac0dffb676de0f4e41e3388ede7d52923f13cd1118"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@cache
def prior_checks():
    if sha(prior.REPORT) != PRIOR_SHA:
        raise ValueError("Pinned S6.4.HR certificate changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return PRIOR_SHA


def residuals():
    groups = {"full_lapse_and_differentiated_Bianchi_identity": flow.identities(),
              "constant_ratio_roots_and_decoupled_exception": flow.root_checks(),
              "sharp_CD_endpoint_mismatch": flow.endpoint_checks(), **controls.checks()}
    for name, values in groups.items():
        if any(sp.simplify(value) != 0 for value in values.values()):
            raise ValueError("A normalized-Hubble identity failed: "+name)
    return groups


def positive_checks():
    bg = flow.bg
    y, c, tau = sp.symbols("y c tau", positive=True)
    inertia = bg.MG2+bg.MF2*y**2
    values = {"inertia": inertia, "dynamic_denominator": 2*inertia**sp.Rational(3, 2),
              "root_denominator": 2*bg.MG2*sp.sqrt(inertia), "second_matter_weight": c*y**3,
              "sharp_CD_endpoint_radius": sp.Rational(8, 5)/tau,
              "interacting_control_NEC_numerator_lower_bound": sp.Rational(2258, 605)}
    if any(value.is_positive is not True for value in values.values()):
        raise ValueError("A stipulated-domain positive expression failed")
    return values


@cache
def build_report():
    previous = prior_checks()
    identities, positive = residuals(), positive_checks()
    arithmetic = fraction_checks.checks()
    control = controls.controls()
    if any(sp.simplify(value) == 0 for value in control.values()):
        raise ValueError("A monotonicity distinction control failed to fire")
    sources = sorted(ROOT.glob("src/p8_bimetric_monotonic/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {"schema": 1, "claim": "P8-S6.5.HR", "date": "2026-09-06",
            "status": "REGULAR_COMMON_FLAT_HR_CONTRACTION_TO_EXPANSION_OBSTRUCTED_INCLUDING_DEGENERATE_TRANSITIONS; GENERAL_PARENT_AND_UV_MATCHING_OPEN",
            "prior_S6_4_HR_sha256": previous,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
            "formulation": "FORMULATION.md", "written_proof": "notes/monotonicity.md",
            "primary_source_audit": "notes/sources.md",
            "primary_sources": {"parent_equations_and_polynomial": "https://arxiv.org/pdf/1111.1655",
                                "interaction_null_anticorrelation": "https://arxiv.org/pdf/1206.3814",
                                "regular_extremality_and_source_hypotheses": "https://arxiv.org/pdf/1211.0214"},
            "exact_residuals": {name: dict.fromkeys(values, "0") for name, values in identities.items()},
            "independent_Fraction_coefficients": arithmetic,
            "positive_expressions": {key: str(value) for key, value in positive.items()},
            "distinction_controls": {key: str(value) for key, value in control.items()},
            "domain": {
                "gravity": "G,F,m4>0; arbitrary constant real beta0,...,beta4; same HR action as pinned S6.4",
                "metric_and_clock": "fixed physical g; positive-root common spatially flat FLRW metrics; Dg=Ng^-1*d/dt",
                "interval": "connected regular interval; smooth positive finite a,b,Ng,Nf at every point; no uniform tail bound or continuation across singular endpoints assumed",
                "matter": "primary g-only NEC matter; extension to independent separately conserved minimally coupled NEC matter on each metric; no shared doubly coupled field",
                "topology": "no restriction on number of root-set components, isolated crossings, branch switches or accumulation points"},
            "nonzero_polynomial_monotone": "P not identically zero: A=G+F*y^2, Z=Hg/sqrt(A), Dg Z<=0 throughout the connected interval",
            "open_dynamic_identity": "-2*A^(3/2)*Dg Z=n_g+c*y^3*n_f on {P(y)!=0}",
            "off_branch_identity": "2*A^(3/2)*Dg Z+n_g+c*y^3*n_f=E_weighted-2F*c*y^3*Df(Hf-Hg/y)",
            "interior_root_identity": "Dg Z=-n_g/(2G*sqrt(A))<=0 where P(y)=0 on an open interval; y is locally constant",
            "boundary_root_argument": "closed root set decomposes into its interior and boundary; every boundary point is approached from {P!=0}; Dg Z is continuous",
            "identically_zero_exception": "P identically zero: Dg Hg=-n_g/(2G)<=0; Z need not be monotone and is not used",
            "sign_conclusion": "for t2>t1, Hg(t1)<=0 implies Hg(t2)<=0; Hg(t1)<0 implies Hg(t2)<0; hence no contraction-to-expansion transition of any degeneracy",
            "sharp_CD_endpoint_mismatch": "max(|Hparent(-tau/2)+8/(5tau)|, |Hparent(tau/2)-8/(5tau)|)>=8/(5tau) in physical g cosmic time; no Hdot error assumption",
            "sharpness": "regular static beta1 Minkowski parent with H=0 attains equality at both endpoints",
            "control_backgrounds": {
                "P_identically_zero": "G=F=m4=1,beta=(3,0,0,0,3),a=e^t,b=e^-t,Ng=Nf=1,no matter: Hg=1,Zdot>0,all four equations exact",
                "P_nonzero_H_not_monotone": "G=F=m4=1,beta=(-2,5/4,0,0,7/4),y=1+t,H=sqrt((7y^2+5/y)/12),c=y+1/H,rho_g=3H^2+2-15y/4,n_g=-rho_g'/(3H),f vacuum; all four equations exact",
                "interacting_control_domain": "9/10<=y<=11/10; n_g>0 and local positive-kinetic canonical reconstruction exists; Hdot(0)=3/8>0 but Zdot(0)=-sqrt(2)/16<0",
                "scope": "exact background controls, not vacuum/perturbation/quantum health claims"},
            "verification_boundary": [
                "Pinned S6.4 and its replayable lineage are checked without edits",
                "Exact symbolic and independent Fraction coefficient identities are replayed; controls bridge to all four pinned Euler equations",
                "Independent covariant/root/control audit is hashed and tested separately",
                "Root-set topology, continuity, mean-value and canonical-reconstruction arguments are written proofs, not formalized PDE or UV results",
                "Primary sources support equations and scope distinctions; no global novelty assertion"],
            "not_established": [
                "Unweighted Hg monotonicity on interacting branches, or Z monotonicity for P identically zero",
                "Healthy arbitrary-beta vacuum or rolling perturbations, a parent matching construction, radiative/interaction control",
                "An exclusion of non-flat or non-common-FLRW/non-bidiagonal metrics, other/singular square-root branches, or singular-endpoint continuation",
                "An exclusion of quantum/NEC-violating/nonminimal/derivative/shared matter, time-dependent interaction coefficients or changed physical frames",
                "A general bimetric or UV-completion exclusion, or completion of S6/P8"]}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Normalized-Hubble certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.5.HR: full root-covered normalized-Hubble monotonicity and sharp CD endpoint mismatch replay passed; general matching OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
