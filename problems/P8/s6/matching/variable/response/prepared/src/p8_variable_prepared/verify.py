"""Read-only analytic prepared-sector and physical source certificate replay."""

import argparse
import hashlib
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

import sympy as sp
from p8_variable_response import verify as prior

from . import analytic, exact, independent, majorants, physical, source

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"analytic-prepared-sector.json"
RESPONSE_SHA = "7ac41e1ca23c9fe9649e6fc703c5d3d236fb599481f921db319416e10865f118"


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
    return (sorted(ROOT.glob("src/p8_variable_prepared/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(prior.REPORT) != RESPONSE_SHA:
        raise ValueError("The frozen full physical-response certificate changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {"S6_21_full_physical_response_and_its_replayed_lineage": RESPONSE_SHA}


def rational_bridge(primary, separate, name):
    if primary.keys() != separate.keys() or any(primary[key] != sp.Rational(separate[key]) for key in primary):
        raise ValueError(f"Independent Fraction bridge failed: {name}")
    return len(primary)


@cache
def independent_bridges():
    count_analytic = rational_bridge(majorants.calibration(), independent.analytic_constants(), "analytic norms")
    count_physical = rational_bridge(physical.calibration(), independent.physical_constants(), "physical norms")
    other = independent.source_constants()
    count_source = rational_bridge(source.coefficient_bounds(), other["coefficients"], "source coefficients")
    switch = source.switch_bounds()
    for name in ("bump_jets", "reciprocal_jets"):
        if switch[name] != [sp.Rational(value) for value in other[name]]:
            raise ValueError("Independent flat-cutoff derivative recurrence failed")
    rational_bridge(switch["zeta_jets"], other["zeta_jets"], "cutoff derivatives")
    if source.source_sup_bound() != sp.Rational(other["source_sup"]):
        raise ValueError("Independent finite source norm failed")
    return {"analytic_constants": count_analytic, "physical_constants": count_physical,
            "source_coefficients": count_source, "flat_cutoff_derivative_orders": 5,
            "all_calibrations_equal": True, "separate_engine_checks": independent.checks()}


def checked_controls():
    delta = sp.Rational(1, 10**12)
    rejected_calls = [lambda value=value: exact.parameters(value) for value in
                      (0, -1, True, 0.000001, sp.Float("0.000001"), sp.Rational(1, 10**8), sp.oo, sp.nan)]
    rejected_calls += [lambda value=value: exact.parameters(delta, value) for value in
                       (-1, 5, True, 1.0, sp.Symbol("K", positive=True))]
    rejected_calls += [lambda: exact.parameters(delta, 0, transfer=True),
                       lambda: analytic.inverse_monomial(0, True),
                       lambda: analytic.inverse_monomial(0, sp.Rational(1, 2)),
                       lambda: physical.fixed_light_error(delta, -1),
                       lambda: physical.fixed_source_error(delta, 1, True),
                       lambda: physical.relative_jet_remainder(delta, sp.Rational(1, 10), "odd"),
                       lambda: physical.relative_jet_remainder(delta, 0, "both"),
                       lambda: source.source_sup_bound(1, {1: 1, 2: 1, 3: 1}),
                       lambda: source.source_sup_bound(1, {True: 1, 2: 1, 3: 1, 4: 1}),
                       lambda: source.source_sup_bound(1, {1: 1, 2: 1, 3: 1, 4: 1.0}),
                       lambda: independent.fraction(True), lambda: independent.fraction(0.1),
                       lambda: independent.inverse_monomial(0, -1), lambda: serialize(0.1)]
    rejected = 0
    for call in rejected_calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(rejected_calls):
        raise ValueError("An inexact, singular or out-of-scope input was admitted")
    d = analytic.derive()
    wrong_inverse = sp.factor(analytic.leading_operator(d["u"]**2/analytic.diagonal(2))-d["u"]**2)
    if wrong_inverse != d["delta"]/88:
        raise ValueError("The omitted delta-shift control changed")
    if physical.fixed_light_error(delta) != 8600*delta or source.source_sup_bound() != 321137485366608000:
        raise ValueError("A quantitative finite-parameter control failed")
    if 3*exact.DELTA_RADIUS/4 <= 300*delta**2:
        raise ValueError("The missing center-vanishing negative control is ineffective")
    return {"rejected_exact_domain_calls": rejected,
            "wrong_unshifted_inverse_residual": "delta/88",
            "K_zero_analytic_branch_admitted": exact.parameters(delta, 0)[1] == 0,
            "K_zero_old_transfer_import_rejected": True,
            "pointwise_sup_bound_alone_does_not_prove_300_delta_squared": True,
            "joint_constant_vanishing_required_for_delta_cubed_symplectic_error": True,
            "locked_linear_coefficient": "4/5", "prepared_linear_coefficient": "14/33",
            "locked_minus_prepared_linear_coefficient": "62/165",
            "center_remainder_at_delta_1e_minus9": serialize(physical.center_remainder(exact.DELTA_MAX)),
            "fixed_light_error_at_delta_1e_minus9": serialize(physical.fixed_light_error(exact.DELTA_MAX)),
            "fourth_cutoff_derivative_retained": True,
            "source_C_sigma_is_finite_not_a_small_stress_assertion": "321137485366608000"}


@cache
def build_report():
    prior_pins = prior_checks()
    p, s = physical.checks(), source.checks()
    residuals = {"exact_full_operator_and_analytic_inverse": analytic.checks(),
                 "physical_center_and_closed_equation": p["residuals"],
                 "pole_cancelled_physical_source": {key: value for key, value in s.items()
                                                   if key != "source_prefactor_margin"}}
    if any(sp.simplify(value) != 0 for group in residuals.values() for value in group.values()):
        raise ValueError("An exact full-operator, inverse, physical or source identity failed")
    if s["source_prefactor_margin"].is_positive is not True:
        raise ValueError("The exact source prefactor bound failed")
    coefficient_margins = majorants.checks()
    audit = ROOT/"tests"/"test_prepared_independent_audit.py"
    if not audit.is_file():
        raise ValueError("The separately authored full-action audit is required before promotion")
    return {
        "schema": 1, "claim": "P8-S6.23.PREPARED", "date": "2026-09-07",
        "status": "CONVERGENT_PREPARED_FULL_TT_SECTOR_AND_PHYSICAL_SOURCE_BOUNDS; ORIGINAL_S6_P8_OPEN",
        "prior_sha256": prior_pins,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
        "formulation": "FORMULATION.md", "written_proof": "notes/proof.md", "source_audit": "notes/sources.md",
        "exact_residuals": {key: dict.fromkeys(group, "0") for key, group in residuals.items()},
        "strict_continuous_margins": serialize({"analytic_coefficients": coefficient_margins,
                                                 "physical_coefficients": p["strict_margins"],
                                                 "source_prefactor": s["source_prefactor_margin"]}),
        "analytic_norm_calibration": serialize(majorants.calibration()),
        "physical_norm_calibration": serialize(physical.calibration()),
        "independent_Fraction_replay": independent_bridges(),
        "checked_controls": checked_controls(),
        "unchanged_physical_parent": {
            "action": "Pinned VARIABLE canonical-clock/free-chi action and exact local background; physical matter metric g",
            "sign_clock": "+---; R_B=-6(DH+2H^2); Einstein action -M^2 R_B/2; u=T/tau; canonical fields divided by M",
            "actual_domain": "M,tau>0; 0<delta<=1e-9; 0<=K=(tau*kcom)^2<=4; |u|<=1/100",
            "separate_transfer_domain": "1<=K<=4, exactly the frozen S6.21 imported uniform transfer band",
            "analytic_device": "Complex |u|<=1/20, |delta|<=1/400, |K|<=4 for the pole-cleared equations only; not actual negative/complex-delta parents",
            "preparation": "Jointly analytic (l,q), Q=Dq, D=(2+delta)(1+u^2)^4-2; even light data(1,0),odd(0,1) at u0",
        },
        "convergent_analytic_sector": {
            "norm": "sum |f_jn(K)| v^j R^n, R=1/20,v=R^2; Banach-valued holomorphy in complex K disk4",
            "inverse": "L0=(delta+8u^2)d_u^2+32u d_u+96; diagonal8(n^2+3n+12); delta shift norm<1/8; inverse<=1/84",
            "derivative_compositions": "L0inv P d_u <=||P||/(70R) for u-degree(P)>=1; second derivative <=||P||/(7R^2) for degree>=2",
            "contraction": "Both exact rational columns are <1/2 in ||l||+v||q||; no formal-only expansion",
            "even_bounds": "||l_e-1||<10v; ||q_e||<3v; q_e(0,0)=0; |q_e(delta,0)|<=3delta",
            "odd_bounds": "||l_o-u||<R/40; ||q_o||<R/5; |q_o,u(delta,0)|<=1/5",
            "leading_q": "q_o=3u/40+weight>=3; q_e=-K(4u^2/55+7delta/2640)+weight>=4",
            "full_Q_remainders": "t=max(sqrt(delta/v),|u|/R): odd <=Dnorm*(R/5)*t^5; even <=3v*Dnorm*t^6",
            "analytic_K_zero": "Constant equal physical g,f is exact at K0; uniqueness gives q_even(K0)=0 and G(K0,u)=0",
            "center_coefficient_l1_tail": "Banach K-Schwarz gives ||q_even/K||<=3v/4; at u0 weight j>=2 gives <=300delta^2 without a geometric-tail factor",
        },
        "physical_prepared_sector": {
            "conserved_symplectic_form": "pi_l=l'-2omega Q,pi_Q=Q'; Omega=1+delta^2*q_e(delta,0)*q_o,u(delta,0); |Omega-1|<=3delta^3/5",
            "physical_Wronskian": "Wg=ge*go'-ge'*go; 3/5<Wg<2 on the full actual slab, including all moving map derivatives",
            "closed_g_equation": "g''+F g'+Gg=0; F=-Wg'/Wg; G=(ge'*go''-ge''*go')/Wg",
            "restricted_positive_kinetic": "9/20<Keff=Omega/Wg<11/6; M^2/(2tau) integral du Keff[(g')^2-Gg^2] represents this source-free subspace at fixed K only",
            "center": "|G(0)/K-1-(14/33)delta|<=48001delta^2, ratio extended analytically to K0",
            "normalization_not_EFT": "Positive restricted symplectic/kinetic normalization is not unrestricted sourced, spatially local EFT health or a cone",
        },
        "fixed_slice_and_actual_source": {
            "inclusion": "Hdelta=Psi0^-1 W^-1 Zanalytic_delta; H0=L=[I2;0]; actual physical Cdelta cancels exactly, not replaced by C0",
            "inclusion_error": "||Hdelta-L||<=200delta; exact full transfer Tdelta Hminus=Hplus",
            "fixed_light": "||Tdelta L-L||<=8600delta from frozen ||Tdelta||<=42; upper rate only, not a sharp leakage coefficient",
            "fixed_old_source": "Output error<=delta[8600||target||+126000000 integral|sigma0|du]; center coefficient is not thereby inherited",
            "retuned_preparation": "f_loaded=zeta F; g_loaded=f_loaded+[(Kf f_loaded')'+Gf K f_loaded]/Utt; full sigma=4/a^3[(Kg g_loaded')'+Gg K g_loaded+Utt(g_loaded-f_loaded)]",
            "pole_cancellation": "P=1/Utt=D/(KR*N) analytic; commutator source uses only zeta derivatives1..4; no inverse-delta source bound",
            "source_support": "Fixed C-infinity flat switch zeta=s(200(u+1/50)-1/2); source inside(-1/50,-1/100), not band-limited",
            "physical_source": "sigma=tau^2 Pi/M^2; probe action+1/2 integral dT a^3 Pi gamma_g; conserved external spatial TT stress",
            "source_coefficients": serialize(source.coefficient_bounds()),
            "genuine_cutoff_calibration": serialize(source.switch_bounds()),
            "source_sup": "C_sigma*(|alpha_even|+|alpha_odd|), C_sigma=321137485366608000; intentionally coarse linear-probe norm",
            "source_L1": "<=C_sigma*(|alpha_even|+|alpha_odd|)/100",
            "retuning": "||sigma_delta-sigma0||sup<=400delta*C_sigma*(|alpha_even|+|alpha_odd|)",
        },
        "not_established": [
            "No regular delta0 full action, change of physical prescription, or original C/D operator/source match",
            "No assumption or theorem that arbitrary incoming states or a fixed generic source select the analytic sector",
            "No sharp or nonzero leading coefficient for fixed-source light-to-heavy leakage",
            "No center14/33 expansion for fixed delta0-source preparation without retuning",
            "No off-shell spatially local EFT for arbitrary low temporal-frequency sources; fixed duration is not a temporal band",
            "No characteristic, causality, positivity or UV verdict from the finite-K center coefficient",
            "No scalar/vector health, quantum state-selection, particle-production or nonlinear source/backreaction theorem",
            "No global background completion and no closure of the original S6/P8 objectives",
        ],
        "verification_boundary": "Written convergent Banach-space proof plus exact full-action identities, independent Fraction arithmetic and separately authored audits; not proof-assistant formalization or numerical mode sampling",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Prepared-sector report differs from exact read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.23.PREPARED: analytic physical/source replay passed; original S6 and P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
