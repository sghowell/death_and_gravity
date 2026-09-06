"""Read-only replay of the scoped composite all-mode principal-cone theorem."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_composite_modes import verify as prior

from . import background, cone, independent

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"composite-cones.json"
PRIOR_SHA = "018295c17af163cc80801e3dc2a9d7309184e4d1f5eca319b2a59ea85f064090"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(value):
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    return str(value)


@cache
def prior_checks():
    if sha(prior.REPORT) != PRIOR_SHA:
        raise ValueError("Pinned S6.7 composite-mode report changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return PRIOR_SHA


@cache
def build_report():
    previous = prior_checks()
    residuals = {"literal_physical_tensor_and_vector_cones": cone.checks(),
                 "full_lapse_background_monotonicity": background.checks()}
    if any(sp.simplify(value) != 0 for values in residuals.values() for value in values.values()):
        raise ValueError("A composite cone or null identity failed")
    rational = independent.checks()
    control = background.coincident_pressure_control()
    if not (control["coincident_lapse_identity"] == control["initial_H"] == 0
            and control["initial_H_prime"] == sp.Rational(1, 6)
            and control["initial_null"] == 1
            and control["initial_Ng"] == control["initial_Nf"] == sp.Rational(1, 2)):
        raise ValueError("The actual zero-vector outside-contract control failed")
    omissions = cone.controls()
    if not (omissions["g_can_exceed_matter_cone"] > 0
            and omissions["f_can_exceed_matter_cone"] > 0
            and omissions["single_metric_beta_zero_does_not_force_c_equal_y"] == 0
            and omissions["omitted_pressure_falsely_keeps_relative_stiffness"] != 0
            and omissions["squared_speed_not_velocity_average"] != 0):
        raise ValueError("A principal-cone omission or domain control failed")
    d = cone.derive()
    sources = sorted(ROOT.glob("src/p8_composite_cones/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1, "claim": "P8-S6.8.COMPOSITE", "date": "2026-09-06",
        "status": "NO_BOUNCE_UNDER_EXPLICIT_ALL_MODE_SUBLUMINAL_TENSOR_AND_STRICT_VECTOR_PRINCIPAL_CONTRACT; GENERAL_EFT_MATCHING_AND_P8_OPEN",
        "prior_S6_7_composite_modes_sha256": previous,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
        "formulation": "FORMULATION.md", "written_proof": "notes/proof.md",
        "source_audit": "notes/sources.md",
        "exact_residuals": {name: dict.fromkeys(values, "0") for name, values in residuals.items()},
        "independent_Fraction_replay": rational,
        "actual_zero_vector_outside_contract_control": encode(control),
        "omission_and_domain_controls": encode(omissions),
        "tensor_principal_cones": {"g_squared": str(d["c_g_squared"]),
                                   "f_squared": str(d["c_f_squared"]),
                                   "positive_velocity_average": "(alpha/r)*cg+(beta*y/r)*cf=1",
                                   "both_subluminal_iff": "c=y"},
        "domain": {
            "geometry": "Smooth regular positive finite common-flat lapses/scales on the positive root, on a connected time interval",
            "action": "Pinned full HR with constant G,F,m4>0, constant real beta_n, constant positive affine alpha,beta; shared canonical scalar",
            "physical_metric": "g_eff=g(alpha I+beta sqrt(g^-1 f))^2,Ne=Ng*s,Ae=a*r,dT=Ne*dt",
            "null_source": "rho+p>=0 on the full source-aware equations, not separate conservation of induced g/f matter",
            "additional_all_mode_formal_principal_contract": "Both tensor speeds<=1 relative to composite matter and regular Xi>0 with vector speed squared>0, at every time",
            "constraint_domain": "k>0 for the relative vector; no Xi=0 kinetic degeneracy or Xi<0 formal high-k ghost included"},
        "implication_chain": [
            "Both tensor cones no faster than composite matter force c=y everywhere",
            "At c=y, mu=y*Q and cVeff^2=mu/Xi; Xi>0 and strict vector speed imply Q>0",
            "Undivided source-aware Q*B=0 implies B=0 everywhere",
            "c=y,B=0 imply ydot=0, constant r, H_eff=H_g/r and Ne=Ng*r",
            "Full g null equation gives dH_eff/dT=-alpha*r*(rho+p)/(2G)<=0",
            "No contraction-to-expansion transition, including arbitrary degenerate zero sets"],
        "pointwise_pressure_branch_consequence": "Q=0 and Xi>0 imply either a tensor cone exceeds matter or the relative vector speed is zero",
        "controls": [
            "Independent full source and physical-clock audits, no main cone-engine import",
            "Actual locally reconstructed Q=0,c=y bounce with hprime=1/6 and zero vector speed fails the strict vector premise",
            "Asymmetric positive-vector local bounce data have one tensor cone outside matter",
            "An exact proportional de Sitter canonical-source solution meets the stated TT/vector subset without bouncing; no scalar health is inferred",
            "Beta=0, one tensor only, omitted source pressure and squared-speed averaging controls are outside or alter the hypotheses"],
        "verification_boundary": "Written characteristic/order/connected-interval/analytic-control proofs plus exact coefficient identities and pinned full lineage; not Lean formalized",
        "not_established": [
            "A general composite, bimetric, DHOST or ultraviolet-completion exclusion",
            "A light-only finite-band EFT verdict or controlled parent matching",
            "A safe rolling heavy threshold or cutoff, quantum front velocity, subcutoff growth, or positivity amplitude",
            "A canonical rolling mass gap identified with mu or vanishing from mu=0 alone",
            "Stability, completeness or physical admissibility of the outside-contract local controls",
            "Completion of S6 or original P8"]}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Composite-cone certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.8.COMPOSITE: scoped all-mode principal-cone no-bounce replay passed; general EFT matching and P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
