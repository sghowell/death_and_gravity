"""Read-only composite-metric background and local-CD certificate replay."""

import argparse
import hashlib
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

import sympy as sp
from p8_bimetric_monotonic import verify as prior

from . import background, bounce, branches, independent, reconstruction

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"composite-background.json"
PRIOR_SHA = "0efd16fa4f02f45f35448e2056fb37c863c29c2ed09203c3351ba9515baff031"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@cache
def prior_checks():
    if sha(prior.REPORT) != PRIOR_SHA:
        raise ValueError("Pinned separate-matter S6.5 certificate changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return PRIOR_SHA


def exact_json(value):
    if isinstance(value, dict):
        return {key: exact_json(item) for key, item in value.items()}
    if isinstance(value, (Fraction, sp.Basic)):
        return str(value)
    return value


def residuals():
    groups = {"literal_lapse_scale_scalar_variation": background.variation_checks(),
              "both_source_aware_undivided_Bianchi_factors": background.bianchi_checks(),
              "branches_and_bounce_intersections": branches.checks(),
              "free_canonical_regular_bounce": bounce.checks(),
              "local_exact_CD_reconstruction": reconstruction.checks(),
              "physical_scale_restoration": reconstruction.scale_checks()}
    for name, values in groups.items():
        if any(sp.simplify(value) != 0 for value in values.values()):
            raise ValueError("A composite-background identity failed: "+name)
    return groups


@cache
def build_report():
    previous = prior_checks()
    groups = residuals()
    rational = independent.checks()
    box = reconstruction.positivity_box()
    positive = branches.positive_factors()
    if any(value.is_positive is not True for value in positive.values()):
        raise ValueError("A stipulated-domain positive branch factor failed")
    controls = bounce.negative_controls()
    if any(sp.simplify(value) == 0 for value in controls.values()):
        raise ValueError("A composite distinction control failed to fire")
    sources = sorted(ROOT.glob("src/p8_composite/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {"schema": 1, "claim": "P8-S6.6.COMPOSITE", "date": "2026-09-06",
            "status": "NEW_PHYSICAL_COMPOSITE_METRIC_HAS_REGULAR_FREE_SCALAR_BOUNCE_AND_LOCAL_EXACT_CD_WITH_RECONSTRUCTED_POTENTIAL; GLOBAL_MATCHING_AND_EFT_HEALTH_OPEN",
            "prior_S6_5_separate_matter_sha256": previous,
            "prior_scope": "Only gravitational conventions, immutable evidence lineage and exclusion boundary are inherited; the separate-matter physical-g theorem is not applied to composite matter",
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
            "formulation": "FORMULATION.md", "written_proof": "notes/background.md",
            "primary_source_audit": "notes/sources.md",
            "primary_sources": {"composite_metric_and_matter_action": "https://arxiv.org/pdf/1408.1678",
                                "generic_constraint_and_EFT_boundary": "https://arxiv.org/pdf/1409.3834",
                                "cosmological_branches": "https://arxiv.org/pdf/1501.02790",
                                "flat_vacuum_cutoff_upper_bound": "https://arxiv.org/pdf/1506.00666",
                                "higher_source_order_trimetric_completion": "https://arxiv.org/pdf/1804.04671"},
            "exact_residuals": {name: dict.fromkeys(values, "0") for name, values in groups.items()},
            "independent_Fraction_coefficients": rational,
            "positive_branch_factors": exact_json(positive),
            "local_CD_exact_rational_box": exact_json(box),
            "nonzero_distinction_controls": exact_json(controls),
            "domain": {
                "action": "HR with G,F,m4>0 and constant real beta_n; one positive-kinetic canonical scalar coupled only to g_eff=g(alpha I+beta sqrt(g^-1 f))^2; alpha,beta>0",
                "physical_frame": "g_eff and its cosmic time T are prescribed physical; neither g nor f alone is the target metric",
                "geometry": "common flat FLRW on the positive square-root branch; smooth positive finite a,b,Ng,Nf",
                "matter": "one shared composite-coupled field; rho+p=(D_eff chi)^2>=0; no separate conservation imposed on the induced g/f sources; no nonnegative-potential assumption"},
            "undivided_Bianchi": "[m4(beta1+2 beta2 y+beta3 y^2)-alpha beta(alpha+beta y)^2 p]*(Ng bdot-Nf adot)=0",
            "dynamic_branch": "On any connected open B=0 branch, Z=Hg/sqrt(G+F y^2) is nonincreasing in g time and H_eff=Hg/(alpha+beta y), excluding contraction-to-expansion there",
            "nondegenerate_bounce_intersections": "At H_eff=B=0 both constituent Hubbles vanish; Q!=0 invokes the open dynamical branch, Q=0 gives explicitly nonpositive effective acceleration. Thus a CD bounce has B!=0 and a pressure-branch neighborhood",
            "pressure_branch": "Q=0 gives Dg Hg=-alpha r^3(rho+p)/(2G), Df Hf=-beta r^3(rho+p)/(2F y^3); individual Hubble monotonicity does not exclude a composite bounce",
            "free_bounce": {
                "parameters": "G=F=M^2,m4=M^2/tau^2,alpha=beta=1,beta_n=(0,0,1,0,0),V=0; formulas use u=T/tau and chi=M varphi",
                "source": "rho_bar=p_bar=2y/(1+y)^2; a_eff=[(1+y)^2/(4y)]^(1/6); a_eff^3 varphi_prime=1",
                "branch": "X=+sqrt(y(2+5y)/3),Y=-sqrt((5+2y)/(3y^2)); D=(y+5)X-(5y+1)Y>0; Ng=-(5y+1)Y/D,Nf=(y+5)X/D,yprime=6y(1+y)XY/D",
                "bounce": "y=1: Ng=Nf=1/2,varphi_prime=1,H_eff=0,dH_eff/dT=7/(36tau^2)>0",
                "double_root": "Q=2y-p_bar(1+y)^2=-(y-1)^2/2 at p_bar=1/2, but yprime=-sqrt(7/3)!=0; no inverse y(p) chart or division by p_y",
                "not_CD_at_any_time_scale": "Normalized a_eff has a_second=7/36,a_fourth=-97/432 and a_fourth-(3/2)*a_second^2=-9/32, whereas every positively time-rescaled normalized CD shape has this invariant zero",
                "scope": "Exact local regular bounce with nonzero free canonical kinetic energy; not CD even after positive time rescaling, not an all-time controlled-curvature or EFT verdict"},
            "local_CD": {
                "parameters": "same gravitational parameters; reconstructed canonical potential, not V=0",
                "target": "a_eff=(1+u^2)^2,Hbar=4u/(1+u^2),u=T/tau; initial y=1,rho_bar=1/2",
                "ODE": "X=sqrt(y^2+(1+y)^3 rho_bar/3),Y=-sqrt(y^-2+(1+y)^3 rho_bar/(3y^3)); yprime=y(1+y)[XY-Hbar(X+Y)]/(X-yY),rho_bar_prime=-3Hbar[rho_bar+2y/(1+y)^2]",
                "lapses": "Ng=[(1+y)Hbar-yY]/(X-yY),Nf=[X-(1+y)Hbar]/(X-yY)",
                "field_and_potential": "varphi_prime=sqrt(rho_bar+p_bar)>0; V(chi(T))=M^2[rho_bar-p_bar]/(2tau^2),chi=M varphi",
                "proven_window": "|T|<=tau/64: 3/4<y<5/4,1/4<rho_bar<3/4,1/7<Ng,Nf<6/7,rho_bar+p_bar>1/2; analytic ODE and strict no-exit bounds",
                "existence_scope": "Constructs one analytic local potential on the traversed field interval and an exact CD solution there; it does not solve a potential prescribed in advance or the full free M1 action"},
            "verification_boundary": [
                "Pinned old lineage replayed without changes; old exclusion hypotheses are explicitly changed",
                "Literal source-aware action variations and independent Fraction coefficients, not only substituted Friedmann guesses",
                "Both equations in each algebraic branch are retained; pressure double roots and bounce intersections never require dividing H, B, Q or p_y",
                "Exact rational bounds plus written analytic ODE continuation/inverse-function proof; no numerical integration used as an existence oracle",
                "Primary sources provide the coupling/branch/cutoff evidence dictionary; no global novelty claim"],
            "not_established": [
                "A global CD solution, globally specified unique reconstructed potential, or full-CD realization by a free M1 scalar",
                "Recovery of the pinned CD/M1 perturbation action in a changed or unchanged physical frame",
                "Full rolling tensor/vector/two-scalar stability, finite-band normalization, interactions or cutoff hierarchy",
                "Generic nonlinear absence of the BD mode or a healthy ultraviolet completion of the bare composite coupling",
                "A computed background-dependent ghost scale; the quoted flat-vacuum upper bound is not a guaranteed safe cutoff here",
                "Smallness of omitted curvature/loop/higher-source terms or inheritance of trimetric completion beyond leading order",
                "S6 matching, positivity closure, quantum cone conclusions, or completion of P8"]}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Composite-background certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.6.COMPOSITE: source-aware branches, free bounce and local exact-CD window replay passed; global/EFT matching OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
