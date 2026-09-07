"""Read-only replay of the literal A26 constant-clock action-domain result."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_uv import verify as prior

from . import coefficients, independent, jets, regularity

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"a26-vacuum-domain.json"
P8 = next(path for path in ROOT.parents if path.name == "P8")
PRIOR_SHA = "6a7ab2a7bf6719d3829795779a876dab0694a5b99383726c14f0880515d2d7d4"
CONTRACT_SHA = "d00df35d5821b05baa61e897725a22ed1ac21d0536b74bdaece3ab0093215901"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@cache
def prior_checks():
    pins = {"S6_1_smooth_extension_context": (prior.REPORT, PRIOR_SHA),
            "adopted_S6_contract": (P8/"s6"/"FORMULATION.md", CONTRACT_SHA)}
    if any(sha(path) != expected for path, expected in pins.values()):
        raise ValueError("An immutable S6 vacuum/context input changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {key: expected for key, (_, expected) in pins.items()}


def bridges():
    data = independent.checks()
    for fixture in data["Laurent_fixtures"]:
        actual = coefficients.residues(sp.Rational(fixture["a"]), sp.Rational(fixture["g"]))
        if [str(actual[key]) for key in ("A3", "A4", "A5")] != fixture["residues"]:
            raise ValueError("Literal symbolic versus Fraction Laurent residue bridge failed")
    bulk = jets.fourier_bulk()
    reconstructed = independent.ibp_bulk()
    if bulk["quartic_after_time_IBP"] != reconstructed["T2_Tp2_after_IBP"]:
        raise ValueError("Direct action versus independent Fourier/log-IBP bridge failed")
    if bulk["fourth_difference_multiplier"] != reconstructed["fourth_difference_multiplier"]:
        raise ValueError("The nonzero fourth finite difference was lost")
    return {"Laurent_fixtures_checked": len(data["Laurent_fixtures"]),
            "direct_quartic_minus_Fraction_log_IBP": "0",
            "fourth_finite_difference_multiplier": "24"}


def omission_controls():
    cancelled = coefficients.residues(sp.Rational(2, 5), sp.Rational(1, 5))
    if cancelled["generic_obstruction"] or any(cancelled[key] for key in ("A3", "A4", "A5")):
        raise ValueError("An isolated r=0 exclusion was promoted to a generic obstruction")
    e = [jets.evaluate_euler((p, 0, 0), (d, d, d)) for p, d in ((1, 1), (2, 1), (3, 2))]
    if e != [2, sp.Rational(1, 2), sp.Rational(16, 9)] or e[2]-e[0]-e[1] != -sp.Rational(13, 18):
        raise ValueError("The Euler additivity countercontrol failed")
    if jets.evaluate_euler((1, 0, 0), (1, 0, 0)) != 0:
        raise ValueError("The homogeneous cancellation control failed")
    if regularity.tilt_data()["physical_instability_established"]:
        raise ValueError("Literal-coordinate differentiability was promoted to instability")
    bad = [lambda: coefficients.residues(True, 0), lambda: coefficients.residues(0.1, 0),
           lambda: coefficients.residues(sp.oo, 0), lambda: coefficients.residues(0, sp.Rational(1, 2)),
           lambda: jets.evaluate_euler((1, 1, 0), (1, 1, 1)),
           lambda: jets.evaluate_euler((1, 2, 0), (1, 1, 1)),
           lambda: jets.evaluate_euler((1, 0), (1, 1, 1)),
           lambda: jets.evaluate_euler((sp.I, 0, 0), (1, 1, 1))]
    rejected = 0
    for call in bad:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(bad):
        raise ValueError("An invalid exact theorem-domain input was admitted")
    return {"r_zero_generic_verdict": False, "degenerate_F2_not_classified": True,
            "homogeneous_Euler": "0", "generic_Euler_outputs": list(map(str, e)),
            "nonadditivity": "-13/18", "A5_pole_small_field_order": 4,
            "printed_A4_regular_difference_not_dropped": True,
            "source_cancellation_checked_before_regularity_claim": True,
            "regularity_is_physical_instability": False, "invalid_inputs_rejected": rejected}


@cache
def build_report():
    pins = prior_checks()
    residuals = {"complete_coefficients_and_analytic_numerator": coefficients.checks(),
                 "full_symmetric_H_Euler_and_nonzero_bulk": jets.checks(),
                 "complete_reconstructed_source_cancellation": regularity.checks()}
    if any(sp.simplify(value) != 0 for group in residuals.values() for value in group.values()):
        raise ValueError("A literal-action, variational or reconstructed-source identity failed")
    files = sorted(ROOT.glob("src/p8_a26_vacuum/*.py"))+sorted(ROOT.glob("tests/*.py"))
    files += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1, "claim": "P8-S6.19.A26", "date": "2026-09-06",
        "status": "LITERAL_A26_GENERIC_CONSTANT_CLOCK_NONQUADRATIC_DOMAIN; MODIFIED_MATCHING_AND_P8_B_OPEN",
        "prior_context_sha256": pins,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in files},
        "formulation": "FORMULATION.md", "written_proof": "notes/proof.md", "source_audit": "notes/sources.md",
        "exact_residuals": {name: dict.fromkeys(group, "0") for name, group in residuals.items()},
        "independent_Fraction_replay": independent.checks(),
        "primary_independent_bridges": bridges(), "checked_omission_controls": omission_controls(),
        "literal_domain": {
            "primary_version": "arXiv:2606.03302v1, 2 June 2026; HTML equations (1)-(4),(17),(19),(23)-(24),(28),(31)",
            "signature": "(-,+,+,+); X=(partial phi)^2 without 1/2; rolling source clock X=-1; constant clock X=0",
            "coefficient_values": "a=a1(phi0),g=g1(phi0),f=1/2-g!=0,r=2g-a!=0 at a defined finite coefficient point",
            "printed_A4": "Literal printed parentheses are default; earlier-Ia comparison differs by 3F2(1-A1^2)/(2(F2-XA1)^2), regular here",
            "restriction": "Fixed Minkowski metric; phi=phi0+epsilon psi with timelike dpsi on a compact local domain",
            "regular_source_points": "Generic phi!=0 points give ordinary lower-derivative Taylor terms; no such expansion is silently assumed for the literal nonsmooth potential at phi=0",
        },
        "theorem": {
            "Laurent_residues_A3_A4_A5": ["r", "-r", "-r^2/(2f)"],
            "small_field_limit": "epsilon^-2 sum_(i=3..5) Ai Li -> r (L3[psi]-L4[psi])/Xpsi; A5 first enters epsilon^4",
            "full_diagonal_Euler": "-2*d0*d1*d2/X, after all six symmetric-H derivatives; embedded 2+1 slice in flat 4D",
            "Euler_counterexample": "E(psi1)=2,E(psi2)=1/2,E(psi1+psi2)=16/9; defect=-13/18",
            "boundary_invariant_identity": "L=div(log|X| J)/2-log|X|*((Boxpsi)^2-Hesspsi^2)/2",
            "bulk_counterexample": "psi=t+l*T(t)*(cos x+cos y); nonzero compact-time T; integrated l^4 coefficient=-int T^2*T'^2/2<0",
            "fourth_difference": "Weights (1,1,-4,-4,6) at k=(2,-2,1,-1,0); quartic multiplier24, identical boundary jets",
            "interpretation": "No bilinear derivative-sector vacuum Hessian for generic r!=0; regular lower-derivative terms cannot cancel it",
            "open_phi_neighborhood": "Z=9(1+z^2)h-2[z^2 tanh(z+.1)+tanh z] is analytic with Z'(0)=-1/20; every open finite interval has regular f!=0,r!=0,phi!=0 points",
            "extension_excluded": "A smooth full action agreeing with the literal action at nonzero X near0 on an open finite-phi neighborhood, even modulo boundary terms",
            "isolated_points": "No claim of the same leading obstruction at r=0; no nondegenerate expansion or healthy Einstein vacuum verdict at f=0",
        },
        "source_regularity": {
            "b": "9/500", "b_choice": "Explicit audit choice, not a printed numerical b; source Eq21 gives n_s=24/25 for epsilon_background=10",
            "absolute_power": "509/250", "Q": "C2 not C3 at phi0=0; second derivative Holder exponent9/250",
            "W2": "C1 not C2 for the stated B'(0)>0 example; W2 itself finite",
            "full_source": "Q[X B'^2-(partial chi)^2]+W2(B^2-chi^2); zero identically for chi=B(phi)",
            "entropy_IBP": "-Q(partial xi)^2-W2 xi^2+2xi[(Q'B'+QB'')(X+1)+QB'(Boxphi+3H)]",
            "background_control": "xi=chi-B(phi); linear source zero at X=-1,Boxphi=-3H; constant B gives W2=0",
            "scope": "Literal-coordinate coefficient differentiability only, not a reduced instability verdict; chi=B(phi) supplies an explicit full-action restricted countercheck at phi0=0,chi0=B(0)",
        },
        "verification_boundary": "Written variational/analytic argument plus exact symbolic, separate Fraction and independently authored audits; not proof-assistant formalized",
        "primary_sources": ["https://arxiv.org/html/2606.03302v1", "https://arxiv.org/html/2501.09985v2"],
        "not_established": [
            "A physical ghost, gradient instability, tachyon or failure of the source's rolling X=-1 perturbative calculation",
            "Absence of every isolated formal stationary point, or the next-order result at an isolated r=0 point",
            "A UV no-completion theorem for a separately named smooth off-tube repair, higher operators or common parent",
            "Equivalence of interacting source chi to frozen canonical free M1",
            "A common-parent V/G/B certificate, quartic scattering reconstruction, finite-gravity Regge/IR or loop bound",
            "Global field-space analyticity of an unknown parent, or a universal DHOST/bounce exclusion",
            "Closure of S6, P8(b), or P8",
        ],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("A26 constant-clock domain certificate differs from exact read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.19.A26: literal generic constant-clock domain replay passed; repaired matching and P8(b) OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
