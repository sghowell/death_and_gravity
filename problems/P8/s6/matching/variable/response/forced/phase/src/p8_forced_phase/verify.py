"""Read-only fixed-pulse phase-separation certificate replay."""

import argparse
import hashlib
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

import sympy as sp
from p8_variable_forced import verify as prior

from . import core, independent

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"fixed-source-phase-separation.json"
FORCED_SHA = "0050320c15694cab2cc40aa8c3e82eaa345283a29d8347c0990e999c669ac755"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serialize(value):
    if isinstance(value, (Fraction, sp.Rational)):
        return str(value)
    if isinstance(value, dict):
        return {key: serialize(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(item) for item in value]
    if value is None or isinstance(value, (str, bool, int)):
        return value
    raise TypeError("Only exact finite report data are admitted")


def source_files():
    return (sorted(ROOT.glob("src/p8_forced_phase/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(prior.REPORT) != FORCED_SHA:
        raise ValueError("The frozen actual-source certificate changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {"S6_26_actual_source_and_recursive_prepared_response_parent_contract": FORCED_SHA}


def controls():
    calls = [lambda value=value: core.calibration(value) for value in
             (True, 0.001, sp.Float(".001"), 0, -1, sp.Rational(1, 100), sp.oo, sp.nan)]
    calls += [lambda value=value: core.calibration(deficit=value) for value in
              (True, 0.0, -1, sp.Rational(1, 99))]
    calls += [lambda value=value: core.separation_bound(momentum_squared=value) for value in
              (True, 0.0, -1, 5)]
    calls += [lambda: core.separation_bound(amplitude=0),
              lambda: core.separation_bound(amplitude=-1),
              lambda: independent.calibration(0.001), lambda: serialize(0.1)]
    count = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            count += 1
    if count != len(calls):
        raise ValueError("An invalid exact-domain input was admitted")
    return {"rejected_inputs": count, "actual_K_zero_admitted": bool(core.separation_bound(momentum_squared=0) > 0),
            "constant_source_class_is_delta_independent": True,
            "C_infinity_sources_exist_with_deficit_at_most_r_over_200": True,
            "actual_subsequence_limits_not_assumed": True,
            "conditional_numeric_bound_does_not_certify_a_profile": True,
            "no_parametrically_low_RMS_proxy_hierarchy_claim": True}


@cache
def build_report():
    pins = prior_checks()
    evidence = core.checks()
    audit = ROOT/"tests"/"test_phase_independent_audit.py"
    if not audit.is_file():
        raise ValueError("The independently authored phase audit is required before promotion")
    bridges = []
    for denominator in (1000, 10000, 10**6):
        for loss in (Fraction(0), Fraction(1, 200), Fraction(1, 100)):
            primary = core.calibration(sp.Rational(1, denominator), sp.Rational(loss))
            separate = independent.calibration(Fraction(1, denominator), loss)
            if primary.keys() != separate.keys() or any(primary[key] != sp.Rational(separate[key]) for key in primary):
                raise ValueError("The separate Fraction error chain disagrees")
            bridges.append({"r": str(Fraction(1, denominator)), "deficit_over_r": str(loss), "comparisons": len(primary)})
    return {
        "schema": 1, "claim": "P8-S6.27.PHASE", "date": "2026-09-07",
        "status": "FIXED_SMOOTH_PHYSICAL_SOURCE_NO_DELTA_RESPONSE_LIMIT; ORIGINAL_S6_P8_OPEN",
        "prior_sha256": pins,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
        "formulation": "FORMULATION.md", "written_proof": "notes/proof.md", "source_audit": "notes/sources.md",
        "exact_residuals": dict.fromkeys(evidence["residuals"], "0"),
        "strict_continuous_margins": serialize(evidence["margins"]),
        "calibration": serialize(core.calibration()),
        "independent_Fraction_bridges": bridges,
        "independent_binomial_and_norm_checks": serialize(independent.checks()),
        "controls": controls(),
        "unchanged_physical_problem": {
            "action_and_source": "Full S6.26/S6.20 coupled TT action; actual g matter frame; sigma=tau^2 Pi/M^2",
            "domain": "Fixed0<r<=1/1000, fixed0<=K<=4;0<delta<=min(1e-9,r^2/100)",
            "data": "All four physical Cauchy data zero at u=-r; measured gamma_g at u=r",
            "pulse": "Fixed delta-independent real0<=sigma<=1, supported[-r,-r/2], integral(1-sigma)<=r/100",
            "smooth_examples": "C-infinity compact pulses inside(-r,-r/2), built with r/800 transition widths",
        },
        "proof_chain": {
            "full_error": "Actual endpoint differs from ag*l0+bg*G0(jH*sigma) by<1200r^4, all coupled/pole terms controlled",
            "light_limit": "The ordinary light Volterra equation has a finite fixed-source delta0 limit",
            "Jost_connection": "c-a-b=-i*mu; zero/one-argument Gauss connection with uniform tails on v in[r/2,r]; |A|^2-|B|^2=1",
            "physical_phase": "mu log(32r^2/delta), mu=sqrt39/2; reflected contribution is delta-independent",
            "amplitude": "C=bg0(r)*sqrt(r)*conj(A)/mu*integral sqrt(v)jH0(v)sigma(-v)exp(i mu log(v/r))dv; |C|>2r^2/25",
            "continuous_moment": "Rectangular moment>5/28; smoothing loss<=1/100; actual source-weight loss<=4r^2; actual map ratio>=1-9r^2",
            "sequences": "delta_n^+=32r^2 exp[-(pi/2-argC+2pi n)/mu], delta_n^-=32r^2 exp[-(3pi/2-argC+2pi n)/mu]",
            "actual_separation": "Finite limsup(g)-liminf(g)>3r^2/20; no assumption that actual g converges on either selected sequence",
        },
        "not_established": [
            "No exclusion of a delta-dependent or memory-retaining effective description",
            "No conclusion for every source, long/approximately bandlimited preparation, selected correlated data or delta-retuned profiles",
            "No small RMS/stiffness-proxy hierarchy, actual rolling spectral gap or cutoff",
            "No nonlinear/quantum/global/scalar/vector, C/D matching or UV verdict",
            "No exact leading asymptotic equality for the actual full response, only a nonzero limsup-liminf lower bound",
            "No original S6 or P8 closure and no frozen ancestor edit",
        ],
        "verification_boundary": "Written uniform special-function/Volterra/limsup proof with exact physical-weight, moment and continuous rational bounds, independent Fraction reconstruction and separately authored audit; not numerical phase sampling or proof-assistant formalization",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The fixed-source phase certificate differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.27.PHASE: fixed-source phase-separation replay passed; original S6/P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
