"""Read-only full physical-TT fixed-slice response certificate replay."""

import argparse
import hashlib
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

import sympy as sp
from p8_variable_beta import verify as prior

from . import bounds, connection, exact, independent, operator, source

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"matched-physical-response.json"
VARIABLE_SHA = "335cd52028baf56b30c75377aa7e6edd7ec86db2799f49b13f467242674fc4d0"
P8 = next(path for path in ROOT.parents if path.name == "P8")
CONTRACT_SHA = "d00df35d5821b05baa61e897725a22ed1ac21d0536b74bdaece3ab0093215901"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serialize(value):
    if isinstance(value, (Fraction, sp.Basic)):
        if isinstance(value, sp.Basic) and value.has(sp.Float, sp.oo, sp.zoo, sp.nan):
            raise TypeError("Only exact finite report expressions are admitted")
        return str(value)
    if isinstance(value, dict):
        return {str(key): serialize(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serialize(item) for item in value]
    if value is None or isinstance(value, (str, int, bool)):
        return value
    raise TypeError("The report must not depend on rounded numeric data")


def source_files():
    files = sorted(ROOT.glob("src/p8_variable_response/*.py"))+sorted(ROOT.glob("tests/*.py"))
    files += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return files


@cache
def prior_checks():
    pins = {"S6_20_actual_variable_background_and_complete_canonical_TT": (prior.REPORT, VARIABLE_SHA),
            "adopted_S6_contract_unchanged": (P8/"s6"/"FORMULATION.md", CONTRACT_SHA)}
    if any(sha(path) != expected for path, expected in pins.values()):
        raise ValueError("A frozen VARIABLE source or adopted matching contract changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {key: expected for key, (_, expected) in pins.items()}


@cache
def independent_bridges():
    primary = bounds.calibration()
    constants = independent.norm_constants()
    if any(primary[key] != sp.Rational(value) for key, value in constants.items()):
        raise ValueError("Independent Fraction versus primary norm reconstruction failed")
    d, src = operator.derive(), source.coefficients()
    expressions = {key: d[key] for key in ("A", "C", "d_cross", "f_cross", "E", "b_analytic")}
    expressions.update({key: src[key] for key in ("j_light", "j_heavy")})
    derivatives = {key: sp.factor(sp.diff(value, d["c"])) for key, value in expressions.items()}
    fixtures = []
    for u0, c0, k0 in ((Fraction(0), Fraction(2), Fraction(1)),
                       (Fraction(-1, 100), Fraction(401, 200), Fraction(2)),
                       (Fraction(1, 10), Fraction(201, 100), Fraction(3, 2))):
        enclosure = independent.profiles(independent.Interval(u0), independent.Interval(c0), independent.Interval(k0))
        subs = {d["u"]: sp.Rational(u0), d["c"]: sp.Rational(c0), d["kbar"]: sp.Rational(k0)}
        for key, expression in expressions.items():
            for index, actual_expr in (((0, 0), expression), ((0, 1), derivatives[key])):
                actual = sp.simplify(actual_expr.subs(subs))
                lower, upper = enclosure[key].get(index).pair()
                if ((actual-sp.Rational(lower)).is_nonnegative is not True
                        or (sp.Rational(upper)-actual).is_nonnegative is not True):
                    raise ValueError("Independent actual value/mixed Taylor derivative enclosure failed")
        fixtures.append({"u": u0, "c": c0, "kbar": k0,
                         "eight_values_and_eight_derivatives_inside_exact_interval_jets": True})
    return {"norm_constants_equal": dict.fromkeys(constants, True),
            "primary_to_independent_actual_coefficient_bridges": serialize(fixtures)}


def checked_controls():
    d = operator.derive()
    center = {d["u"]: 0, d["c"]: 2}
    if sp.simplify(d["C"].subs(center)) != sp.Rational(48, 5) or d["omega"].subs(center) != 0:
        raise ValueError("The omega-prime omission control changed")
    if bounds.transfer_error(exact.DELTA_MAX) != 40:
        raise ValueError("The large-box error must not be misreported as small")
    small = sp.Rational(1, 10**21)
    if bounds.transfer_error(small) != sp.Rational(1, 250):
        raise ValueError("The explicit nontrivial finite-delta example failed")
    bad = [lambda value=value: exact.parameters(value) for value in
           (0, -1, True, 0.000001, sp.Float("0.000001"), sp.Rational(1, 10**8), sp.oo, sp.nan)]
    bad += [lambda value=value: exact.parameters(small, value) for value in
            (0, 3, True, 1.0, sp.Symbol("k", positive=True))]
    bad += [lambda: bounds.source_error(small, -1), lambda: bounds.prepared_error(small, True, 1),
            lambda: bounds.outer_picard_tail(sp.Rational(3, 100), 1),
            lambda: bounds.outer_picard_tail(0, True), lambda: bounds.exponential_majorant(1),
            lambda: independent.Interval(True), lambda: independent.Interval(0.1),
            lambda: independent.Interval(-1, 1).inverse(),
            lambda: independent.Interval(-1).sqrt(), lambda: serialize(0.1)]
    rejected = 0
    for call in bad:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(bad):
        raise ValueError("An inexact, singular or out-of-scope domain input was admitted")
    return {"all_invalid_exact_domain_calls_rejected": rejected,
            "omega_zero_does_not_remove_C_center": "48/5",
            "mass_remainder_center_limit": "-32", "mass_remainder_punctured_limit": "-141",
            "full_mass_remainder_joint_analyticity": False,
            "whole_box_error_upper_is_not_small": "40",
            "delta_1e_minus21_transfer_error_upper": "1/250",
            "physical_center_heavy_to_light_source_ratio": "-2",
            "prepared_light_error_retains_source_L1": True,
            "B_zero_still_leaves_transmission_phase": True,
            "power_basis_determinant_minus_one_is_not_Cauchy_symplectic_test": True}


@cache
def build_report():
    previous = prior_checks()
    coefficient = bounds.coefficient_checks()
    residuals = {"full_operator_and_both_chart_maps": operator.checks(),
                 "hypergeometric_connection_and_phase": connection.checks(),
                 "actual_g_source_and_g_jet_observability": source.checks(),
                 "mass_pole_and_nonanalytic_path_controls": coefficient["residuals"]}
    if any(sp.simplify(value) != 0 for group in residuals.values() for value in group.values()):
        raise ValueError("An exact operator, chart, connection or physical-source identity failed")
    margins = bounds.checks()
    own = independent.checks()
    audit = ROOT/"tests"/"test_response_independent_audit.py"
    if not audit.is_file():
        raise ValueError("The separately authored physical-action audit must exist before freezing")
    return {
        "schema": 1, "claim": "P8-S6.21.RESPONSE", "date": "2026-09-06",
        "status": "FULL_FIXED_SLICE_PHYSICAL_TT_TRANSFER_AND_SOURCE_PREPARATION_BOUNDS; ORIGINAL_S6_P8_OPEN",
        "prior_sha256": previous,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
        "formulation": "FORMULATION.md", "written_proof": "notes/proof.md", "source_audit": "notes/sources.md",
        "exact_residuals": {key: dict.fromkeys(group, "0") for key, group in residuals.items()},
        "continuous_coefficient_margins": serialize(coefficient["strict_continuous_margins"]),
        "norm_calibration": serialize(margins), "independent_Fraction_interval_replay": serialize(own),
        "independent_primary_bridges": independent_bridges(), "checked_controls": checked_controls(),
        "unchanged_action_and_physical_frame": {
            "parent": "Exact S6.20 variable-beta canonical-clock/free-chi action and actual local solution; not an original C/D operator match",
            "signature_curvature": "+---; R_B=-6(DH+2H²); Einstein action -M²R_B/2; actual physical metric g",
            "clock_and_coefficients": "u=T/tau, c=2+delta, a=(1+u²)²,b=2/(1+u²)²; each finite positive delta is an actual member of the pinned family",
            "domain": "M,tau>0; 0<delta<=1e-9, 1<=kbar=tau*kcom<=2; slices u=+-1/100; optional source support(-2/100,-1/100)",
            "singular_limit": "delta=0 only defines regular punctured outer operators; it is not a regular full parent at u=0",
        },
        "exact_operator_and_norm": {
            "full_equations": "l''+A*l+C*Q+u*d_cross*Q'=j_l; Q''+B*Q+E*l+u*f_cross*l'=j_Q, including both moving-weight adjoints",
            "regularity_dictionary": "b_analytic=B-mass is analytic; b_delta=B-80/(delta+8u²) is only bounded, with different mass-only path limits -32,-141",
            "state_norm": "Euclidean (l,l',Q/sqrt(r),(uQ'-Q/2)/(mu sqrt(r))); induced operator norm; mu=sqrt(39)/2",
            "actual_endpoint_map": "P_delta^side=C_delta(side*a*) W_side(a*) Psi_side(a*); no delta0 physical-map substitution",
            "physical_clock": "C_delta maps canonical fields divided by M to (gamma_g,tau*gamma_g,T,gamma_f,tau*gamma_f,T); raw derivatives keep explicit tau factors",
            "outer_waves": "Full delta0 coupled punctured Volterra solution, ||R||<=64r, both wave norms<=exp(64r); factorial Picard tail exposed",
        },
        "matched_transfer_theorem": {
            "charts": "m=delta^(1/3), epsilon=sqrt(delta), eta=delta/(8m²)=m/8, z=asinh(sqrt8*u/epsilon); overlapping at u=+-m",
            "inner_exact": "Q=sqrt(r_e)*psi,r_e=sqrt(u²+delta/8); psi_zz+[39/4+3/4 sech²z]psi=0 is the universal block; full coupled remainder norm<=72r_e",
            "inner_L1": "integral||remainder|| dz<=144m; full propagator<4; direct Duhamel error<=4*2*144m=1152m",
            "outer_difference": "||M_delta-M0||<=2delta/r², all mass variation treated separately; each normalized outer transfer differs from I by<=8m",
            "boundary_and_tails": "Exact fractional-weight overlap differs from sign map by<=m/8; both Jost tails and both physical phase errors together<=2m",
            "operator_error": "||(P_delta^+)^-1 T_phys_delta P_delta^- -diag(I2,S_epsilon_real)||<=40000delta^(1/3)",
            "finite_example": "At delta<=1e-21 this error is <=1/250; at delta=1e-9 the bound40 is not small",
            "endpoint_norms": "||P_delta^side||<3000, ||(P_delta^side)^-1||<280, in the dimensionless actual physical Cauchy clock",
        },
        "exact_connection": {
            "A": serialize(connection.coefficients()["A"]), "B": "i/sinh(pi*mu)",
            "S_power": "[[-conj(B),conj(A)*(4sqrt2/epsilon)^(2i mu)],[A*(epsilon/(4sqrt2))^(2i mu),-B]]",
            "normalization": "Unitary complex power coordinates Z+-=(C1-+i C2)/sqrt2; real data have conjugate pair; left positive-frequency exponent is opposite right",
            "moduli": "|A|²-|B|²=1; 0<|B|<1/1000; ||S||<2; power-coordinate determinant=-1 from radial orientation",
            "phase_control": "A carries a nonvanishing log-epsilon transmission phase even if B is set to zero; no reflection-only or quantum-production inference",
        },
        "physical_source_and_data": {
            "conserved_probe": "delta T^ij=a^-2 Pi(u)e^ij cos(kcom z), delta T^0mu=0, unit transverse traceless polarization; arbitrary smooth time Pi",
            "literal_probe_action": "+1/2 integral dT a³ Pi gamma_g; sigma=tau²Pi/M²; normalized forcing a³sigma/(2fsum),-a³w2sigma/(2fR)",
            "loading_operator": "gamma_f=zeta gamma_f,target; gamma_g=gamma_f+[(Kf gamma_f')'+Gf kbar² gamma_f]/U_TT; sigma0=4/a³[(Kg gamma_g')'+Gg kbar² gamma_g+U_TT(gamma_g-gamma_f)]",
            "fixed_source": "Full delta0 punctured homogeneous target and a fixed smooth flat cutoff on(-2a*,-a*); one fixed sigma0 for every positive delta; retarded zero initial data",
            "source_error": "Incoming outer coefficient error<=3e6 delta S, S=integral|sigma0|du, finite and explicit; no invented S value",
            "paired_targets": "C=(0,0,1,0) and C=(1,0,0,0), relative and prepared regular-light outer data, respectively",
            "outgoing_error": "40000delta^(1/3)||C||+(2+40000delta^(1/3))*3e6 delta S",
            "prepared_light": "Vanishing bounded relative leakage, not a sharp O(delta) or certified low-temporal-band result",
            "g_only_observability": "After the source ends, (g,g',f,f')->(g,g',g'',g''') has determinant (U_TT/Kg)²>0 at fixed punctured slices; raw proper jets retain tau^-j",
            "relative_response": "For nonzero prepared relative data, two log-phase subsequences give different physical outputs; finite endpoint maps cannot erase full/g-four-jet nonconvergence",
        },
        "not_established": [
            "No delta0 regular full action or uniform global-background continuation",
            "No reduced scalar/vector health, temporal spectral gap or heavy-state suppression supplied by the physical model",
            "No sharp prepared-light leakage order or controlled omission error below the O(delta) center locked-cone excess",
            "Fixed source duration and spatial k band are not a low temporal band; a compact pulse is not band-limited",
            "No nonlinear finite-amplitude source/backreaction, vacuum-production or particle-number claim",
            "No original DHOST C/D operator/source match, light-only EFT exclusion, positivity or UV-completion verdict",
            "No change to completed scoped photon objective, adopted S6 contract or original P8 open status",
        ],
        "verification_boundary": "Exact action/chart/connection identities, full continuous Fraction interval jets, independent norm arithmetic and separately authored audits support the written linear-ODE estimates; no proof-assistant formalization or numerical propagator sampling substituted for the proof",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Matched physical response certificate differs from read-only exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.21.RESPONSE: full fixed-slice physical TT/source replay passed; original S6 and P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
