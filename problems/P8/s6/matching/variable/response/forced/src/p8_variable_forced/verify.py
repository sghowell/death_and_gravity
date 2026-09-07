"""Read-only replay of the actual-source causal nonlocal response gate."""

import argparse
import hashlib
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

import sympy as sp
from p8_variable_prepared import verify as prior

from . import bounds, coefficients, exact, independent, kernel, operator, source

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"causal-physical-response.json"
PREPARED_SHA = "a8d89126980dd8ec458b3161bc8b397f74211736f02e64d4e9e29b95e43a0ade"
CONTRACT_SHA = "d00df35d5821b05baa61e897725a22ed1ac21d0536b74bdaece3ab0093215901"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serialize(value):
    if isinstance(value, (Fraction, sp.Basic)):
        if isinstance(value, sp.Basic) and value.has(sp.Float, sp.oo, sp.zoo, sp.nan):
            raise TypeError("Report arithmetic must be exact and finite")
        return str(value)
    if isinstance(value, dict):
        return {str(key): serialize(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(item) for item in value]
    if value is None or isinstance(value, (str, int, bool)):
        return value
    raise TypeError("Unsupported or inexact report data")


def source_files():
    return (sorted(ROOT.glob("src/p8_variable_forced/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(prior.REPORT) != PREPARED_SHA:
        raise ValueError("The frozen prepared-sector certificate changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    contract = next(path for path in ROOT.parents if path.name == "s6")/"FORMULATION.md"
    if sha(contract) != CONTRACT_SHA:
        raise ValueError("The adopted original S6 contract changed")
    return {"S6_23_prepared_and_replayed_S6_21_S6_20_lineage": PREPARED_SHA,
            "adopted_original_S6_contract": CONTRACT_SHA}


def bridge(primary, separate, name):
    if primary.keys() != separate.keys() or any(primary[key] != sp.Rational(separate[key]) for key in primary):
        raise ValueError(f"The independent Fraction bridge failed: {name}")
    return len(primary)


@cache
def independent_checks():
    rows = []
    for denominator in (100, 1000, 10**6):
        rr = sp.Rational(1, denominator)
        ff = Fraction(1, denominator)
        rows.append({"r": str(rr),
                     "coefficient_constants": bridge(coefficients.envelope(rr), independent.coefficient_constants(ff), "coefficients"),
                     "causal_constants": bridge(bounds.calibration(rr), independent.response_constants(ff), "response")})
    # The arbitrary-r theorem uses the coefficientwise positivity replay,
    # not agreement at these three arithmetic regression fixtures alone.
    return {"exact_arithmetic_regression_rows": rows, "continuous_polynomial_replay": independent.checks()}


def domain_controls():
    delta = sp.Rational(1, 10**12)
    calls = [lambda value=value: exact.parameters(value) for value in
             (True, 0.0, sp.Float("1e-12"), 0, -1, sp.Rational(1, 10**8), sp.oo, sp.nan)]
    calls += [lambda value=value: exact.parameters(delta, momentum_squared=value) for value in
              (True, 0.0, -1, 5, sp.Symbol("K", positive=True), sp.zoo)]
    calls += [lambda value=value: exact.radius(value) for value in
              (True, 0.01, 0, -1, sp.Rational(1, 50), sp.oo)]
    calls += [lambda: exact.parameters(sp.Rational(1, 10**9), sp.Rational(1, 10**6)),
              lambda: bounds.physical_response_error(delta, -1),
              lambda: kernel.universal_kernel(sp.Rational(1, 50), 0, delta),
              lambda: independent.rational(True), lambda: independent.rational(0.1),
              lambda: serialize(0.1)]
    count = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            count += 1
    if count != len(calls):
        raise ValueError("An invalid or out-of-domain input was admitted")
    if bounds.physical_response_error(delta, 1, momentum_squared=0) != sp.Rational(1, 10**5):
        raise ValueError("The actual K-zero source error calibration failed")
    if bounds.physical_response_error(delta, 1) != bounds.physical_response_error(delta**2, 1):
        raise ValueError("A false delta-convergent interpretation entered the bound")
    d = operator.derive()
    if sp.factor(d["ag"]*d["jL"]-2) == 0:
        raise ValueError("The omitted heavy source must not reproduce physical normalization")
    return {"rejected_exact_domain_calls": count, "actual_K_zero_admitted": True,
            "retarded_kernel_zero_for_output_before_source": kernel.universal_kernel(-sp.Rational(1, 200), 0, delta) == 0,
            "mass_remainder_center_path": "-32", "mass_remainder_punctured_path": "-141",
            "b_delta_not_jointly_continuous_or_uniform_delta_differentiable": True,
            "omitted_jH_source_fails_actual_g_diagonal": True,
            "actual_full_source_diagonal": "2", "prepared_limit_source_diagonal": "2/5",
            "complementary_diagonal_limit": "8/5",
            "impulse_diagonal_not_low_frequency_error_verdict": True,
            "fixed_r_error_bound_does_not_decrease_with_delta": True,
            "RMS_over_stiffness_proxy_lower": "6/13",
            "RMS_is_not_an_exact_bandlimit_and_proxy_is_not_a_rolling_gap": True}


@cache
def build_report():
    pins = prior_checks()
    cc = coefficients.checks()
    groups = {"literal_action_physical_source": operator.checks(),
              "positive_mass_pole_isolation": cc["residuals"],
              "actual_clock_retarded_kernel": kernel.checks(),
              "symplectic_source_projection": source.checks()}
    if any(sp.simplify(value) != 0 for group in groups.values() for value in group.values()):
        raise ValueError("An exact full-action, clock or source residual failed")
    audit = ROOT/"tests"/"test_forced_independent_audit.py"
    if not audit.is_file():
        raise ValueError("The separately authored source-aware audit is required before promotion")
    return {
        "schema": 1, "claim": "P8-S6.26.FORCED", "date": "2026-09-07",
        "status": "CAUSAL_NONLOCAL_PHYSICAL_G_RESPONSE_BOUND; ORIGINAL_S6_P8_OPEN",
        "prior_sha256": pins,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
        "formulation": "FORMULATION.md", "written_proof": "notes/proof.md", "source_audit": "notes/sources.md",
        "exact_residuals": {name: dict.fromkeys(group, "0") for name, group in groups.items()},
        "strict_continuous_margins": serialize({"literal_coefficients": cc["strict_margins"],
                                                 "causal_source_norms": bounds.checks()}),
        "coefficient_envelope_at_r_1_over_100": serialize(coefficients.envelope()),
        "causal_calibration_at_r_1_over_100": serialize(bounds.calibration()),
        "independent_Fraction_replay": serialize(independent_checks()),
        "checked_controls": domain_controls(),
        "unchanged_physical_parent": {
            "action": "Pinned VARIABLE full two-TT action, not a locked ansatz; source/matter metric g",
            "clock_convention": "+---; R_B=-6(DH+2H^2); Einstein -M^2 R_B/2; u=T/tau; canonical fields divided by M",
            "domain": "M,tau>0; 0<r<=1/100; I=[-r,r]; 0<delta<=min(1e-9,r^2/100); 0<=K=(tau*kcom)^2<=4",
            "initial_data": "All four physical/canonical Cauchy data zero at -r; no implicit analytic or quantum state selection",
            "actual_source": "sigma=tau^2 Pi/M^2; +1/2 integral dT a^3 Pi gamma_g; conserved external spatial TT polarization",
            "source_class": "Arbitrary real L-infinity(I), norm S; smooth flat-past pulses included; fixed spatial K per mode",
        },
        "causal_reduction": {
            "exact_schur": "(L_L-Bopstar GH Bop)l=jL sigma-Bopstar GH(jH sigma); Q=GH(jH sigma-Bop l)",
            "inverses": "All retarded with zero data at -r; GL=(d_u^2+A)^-1; GH=(d_u^2+B)^-1",
            "universal_kernel": "G0=sqrt(rho(u)rho(s))*Im[f+(z(u))*conj f+(z(s))]/mu for u>=s, zero otherwise",
            "finite_clock": "rho=sqrt(u^2+delta/8), z=asinh(sqrt8*u/sqrt(delta)), mu=sqrt39/2",
            "Jost": "f+(z)=exp(i mu z)*2F1(-1/2,3/2;1-i mu;(1-tanh z)/2); W=-2i mu; positive retarded jump1",
            "heavy_norms": "GH:Linf->C0 <=101r^2/75; Bopstar GH:Linf->Linf <64r^2; GH-G0:Linf->C0 <292r^4",
            "light_norm": "N(l)=||l||inf+r||l'||inf; ||GL||Linf->N<5r^2; ||Bop||N->Linf<=10+60r^2",
            "feedback": "<=320r^4(10+60r^2)<1/10000; exact causal Neumann series",
            "approximation": "l0=GL(jL sigma); Qapp=G0[jH sigma-Bop l0]; gapp=ag*l0+bg*Qapp with actual delta-dependent map",
            "errors": "N(l-l0)<323r^4S; ||Q-Qapp||inf<296r^4S; ||g-gapp||inf<915r^4S<=1000r^4S",
            "fixed_example": "r=1/100 gives <=1e-5*S for every admitted delta, including K0; absolute C0 source-to-field bound only",
        },
        "source_and_scale_boundaries": {
            "prepared_projection": "Exact symplectic projection has closed-g source a^3 Wg/(2Omega)=a^3/(2Keff), not full source2",
            "center_limit": "Projected retarded diagonal tends2/5, full diagonal2, complement8/5; not low-band error lower bound",
            "RMS_class": "Nonzero sigma in H1_0(I); Omega_sigma=||sigma'||2/||sigma||2>=pi/(2r)",
            "stiffness_proxy": "m_proxy=inf sqrt(B)>0; r^2 B>=8000/801-161r^2>9; r^2 B(r)<(13/4)^2",
            "declared_hierarchy": "Omega_sigma/m_proxy>6/13; no parametrically small version of this particular RMS/proxy hierarchy via delta",
            "state_dependence": "Nonzero initial heavy data add a homogeneous causal response; not silently projected out",
        },
        "not_established": [
            "No local-in-time derivative expansion or arbitrary-source local light-only EFT",
            "No delta-convergence of the stated fixed-r error bound, or relative error where the response may vanish",
            "No error budget resolving an O(delta) locked-cone effect or importing the prepared center coefficient for arbitrary forcing",
            "No uniform raw physical derivative error from the C0 norm or a low-frequency error from the impulse diagonal",
            "No identification of the declared RMS scale with a Fourier bandlimit or sqrt(B) with a rolling spectral gap/cutoff",
            "No exclusion of other source/preparation/history classes or every possible low-energy reduction",
            "No characteristic, causality-violation, positivity, UV, scalar/vector health or original C/D matching verdict",
            "No nonlinear source/backreaction, quantum state-selection, particle-production or global background theorem",
            "No change to the frozen physical prescription and no original S6/P8 closure",
        ],
        "verification_boundary": "Full written retarded-ODE proof, exact literal source/clock identities, continuous positive-polynomial envelopes, independent Fraction replay and separately authored audits; not numerical mode sampling or proof-assistant formalization",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The forced physical-response report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.26.FORCED: causal nonlocal actual-source replay passed; original S6 and P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
